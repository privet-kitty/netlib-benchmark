(eval-when (:compile-toplevel :load-toplevel :execute)
  (ql:quickload '(:cp/lp :cp/sparse-simplex :cl-mps)))

(defpackage :netlib-benchmark
  (:use :cl :cp/lp :cp/sparse-simplex)
  (:import-from :cl-mps)
  (:import-from :cp/csc #:csc-float)
  (:local-nicknames (#:m #:cl-mps)))
(in-package :netlib-benchmark)

(setq *default-pathname-defaults*
      (uiop:pathname-directory-pathname (uiop:current-lisp-file-pathname)))
(uiop:chdir *default-pathname-defaults*)
(defparameter *netlib-dir* (merge-pathnames #P"netlib/" *default-pathname-defaults*))

(defun make-master-table (paths)
  (let ((table (make-hash-table :test #'equal)))
    (dolist (path paths)
      (uiop:with-input-file (in path)
        (read-line in)
        (loop for line = (read-line in nil nil)
              while line
              for (name nz size obj) = (uiop:split-string line :separator '(#\,))
              do (setf (gethash name table)
                       (let ((*read-default-float-format*
                               (sb-ext:typexpand 'csc-float)))
                         (read-from-string obj))))))
    table))

(defparameter *table* (make-master-table (directory (merge-pathnames  "table-netlib.csv" *netlib-dir*))))

(defparameter *timeout* 1000d0)

(defun convert-problem (problem)
  (let ((res (make-lp-problem :sense (m:problem-sense problem)))
        (mps-vars (m:problem-variables problem))
        (mps-constrs (m:problem-constraints problem))
        (mps-obj (m:problem-objective problem))
        (vars (make-hash-table :test #'equal)))
    ;; vars
    (maphash (lambda (name mps-var)
               (let ((var (new-lp-var res
                                      (m:var-lo mps-var)
                                      (m:var-up mps-var)
                                      (m:var-name mps-var))))
                 (setf (gethash name vars) var)))
             mps-vars)
    ;; constraints
    (loop for mps-constr being each hash-value of mps-constrs
          for expr = (linear-expr)
          for lo = (m:constraint-lo mps-constr)
          for up = (m:constraint-up mps-constr)
          do (maphash (lambda (var-name coef)
                        (let ((var (gethash var-name vars)))
                          (linear-expr-add-var expr coef var)))
                      (m:constraint-coefs mps-constr))
             (new-lp-constr res expr lo up))
    ;; objective
    (let ((obj (linear-expr))
          (obj-offset (- (coerce (m:constraint-rhs mps-obj) 'double-float))))
      (maphash (lambda (var-name coef)
                 (let ((var (gethash var-name vars)))
                   (linear-expr-add-var obj coef var)))
               (m:constraint-coefs mps-obj))
      (setf (lp-problem-objective res) obj)
      (values res obj-offset))))

(defun proc-instance (path)
  (with-open-file (in path)
    (let* ((mps-problem (let* ((*read-default-float-format* (sb-ext:typexpand 'csc-float)))
                          (m:read-fixed-mps in)))
           (name (string-trim '(#\Space) (m:problem-name mps-problem)))
           (master-obj-value (gethash name *table*)))
      (m:problem-name mps-problem)
      (unless master-obj-value
        (return-from proc-instance
          (values nil nil nil nil nil)))
      (multiple-value-bind (problem obj-offset) (convert-problem mps-problem)
        (let* ((start-time (get-internal-real-time))
               (status (handler-bind
                           ((sb-ext:timeout
                              (lambda (c) (declare (ignorable c))
                                (return-from proc-instance
                                  (values name :not-solved nil nil "inf"))))
                            (error
                              (lambda (c) (declare (ignorable c))
                                (return-from proc-instance
                                  (values name :not-solved nil nil nil)))))
                         (sb-ext:with-timeout *timeout*
                           (lp-problem-solve problem #'slp-dual-primal!))))
               (end-time (get-internal-real-time))
               (obj-value (+ obj-offset (lp-problem-obj-value problem)))
               (time (* (- end-time start-time) 1d-3)))
          (values name
                  status
                  (float obj-value 1d0)
                  (float master-obj-value 1d0)
                  time))))))

(defun proc-all (&optional (stream *standard-output*))
  (format stream "name,status,obj,master_obj,time~%")
  (dolist (path (directory (merge-pathnames "*.SIF" *netlib-dir*)))
    (print path)
    (multiple-value-bind (name status obj master-obj time) (proc-instance path)
      (when name
        (let ((*read-default-float-format* 'double-float))
          (format stream "~A,~A,~A,~A,~A~%" name (string status) obj master-obj time)
          (finish-output stream))))))

#-swank
(uiop:with-output-file (out "result.csv" :if-exists :supersede)
  (proc-all out))
