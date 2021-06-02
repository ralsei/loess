#lang typed/racket
(require math/array math/flonum math/matrix
         (only-in pfds/deque/bankers
                  deque deque->list enqueue enqueue-front Deque head+tail)
         (prefix-in d: pfds/deque/bankers))

(: matrix-min* (-> (Matrix Flonum) Flonum))
(define (matrix-min* mat)
  (apply min (matrix->list mat)))

(: matrix-max* (-> (Matrix Flonum) Flonum))
(define (matrix-max* mat)
  (apply max (matrix->list mat)))

(: matrix-argmin (-> (Matrix Flonum) Indexes))
(define (matrix-argmin mat)
  (define-values (m n) (matrix-shape mat))
  (define indices (indexes-array (vector m n)))
  (for/fold : Indexes
      ([current-min-value : Flonum +inf.0]
       [current-min-index : Indexes #(0 0)]
       #:result current-min-index)
      ([val (in-array mat)]
       [idx (in-array indices)])
    (if (< val current-min-value)
        (values val idx)
        (values current-min-value current-min-index))))

(: matrix-argmax (-> (Matrix Flonum) Indexes))
(define (matrix-argmax mat)
  (define-values (m n) (matrix-shape mat))
  (define indices (indexes-array (vector m n)))
  (for/fold : Indexes
      ([current-max-value : Flonum -inf.0]
       [current-max-index : Indexes #(0 0)]
       #:result current-max-index)
      ([val (in-array mat)]
       [idx (in-array indices)])
    (if (> val current-max-value)
        (values val idx)
        (values current-max-value current-max-index))))

(: tricubic (-> Flonum Flonum))
(define (tricubic n)
  (cond [(> (abs n) 1) 0.0]
        [else (cast (expt (expt (- 1.0 (abs n)) 3.0) 3.0) Flonum)]))

(: get-indexer (-> (Matrix Flonum) Integer (U Slice (Sequenceof Integer))))
(define (get-indexer distances window)
  (match-define (vector _ idx) (matrix-argmin distances))
  (define-values (_ n) (matrix-shape distances))

  ; TODO: this could be wildly more efficient
  (cond [(zero? idx) (:: 0 window)]
        [(= idx (sub1 n)) (:: (- n window) n)]
        [else (for/fold : (Sequenceof Integer)
                  ([deq : (Deque Integer) (deque idx)] #:result (sort (deque->list deq) <))
                  ([_ (in-range window)])
                (define h (d:head deq))
                (define t (d:last deq))
                (cond [(zero? h) (enqueue (add1 t) deq)]
                      [(= t (sub1 n)) (enqueue-front (sub1 t) deq)]
                      [(< (matrix-ref distances 0 (sub1 h))
                          (matrix-ref distances 0 (add1 t)))
                       (enqueue-front (sub1 h) deq)]
                      [else (enqueue (add1 t) deq)]))]))

(: get-weights (-> (Matrix Flonum) (U Slice (Sequenceof Integer)) (Matrix Flonum)))
(define (get-weights distances indexer)
  (define indexed (array-slice-ref distances (list (::) indexer)))
  (define max-dist (matrix-max* indexed))
  (array-map (λ ([v : Flonum]) (tricubic (/ v max-dist))) indexed))

(: normalize-x (-> (Matrix Flonum) Flonum Flonum))
(define (normalize-x xs x)
  (define xs-min (matrix-min* xs))
  (define xs-max (matrix-max* xs))
  (/ (- x xs-min) (- xs-max xs-min)))

(: denormalize-y (-> (Matrix Flonum) Flonum Flonum))
(define (denormalize-y ys y)
  (define ys-min (matrix-min* ys))
  (define ys-max (matrix-max* ys))
  (+ (* y (- ys-max ys-min)) ys-min))

(: loess-fit (->* [(Vectorof Flonum) (Vectorof Flonum) Integer]
                  [#:degree Positive-Integer]
                  (-> Flonum Flonum)))
(define (loess-fit xs ys window #:degree [degree 1])
  (define xs-mat (matrix-normalize (->row-matrix xs)))
  (define Y (matrix-normalize (->row-matrix ys)))

  (λ (x)
    (define norm-x (normalize-x xs-mat x))
    (define distances (array-map (λ ([v : Real]) (abs (- v norm-x))) xs-mat))
    (define indexer (get-indexer distances window))
    (define weights (get-weights distances indexer))

    ; β = (X^T*WX)^(-1)(X^T*WY), linearly
    (define W (matrix* weights (identity-matrix window)))

    (define X^n (matrix-expt X degree))
    (define X^T (matrix-transpose X^n))
    (define X^T*W (matrix* X^T W))
    (define β (matrix* (matrix-inverse (matrix* X^T*W X^n))
                       (matrix* X^T*W Y)))
    (pretty-print β)

    3.0))
