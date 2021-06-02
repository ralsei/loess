#lang typed/racket
(require math/array math/matrix
         plot/utils
         plot)

(: vmin (-> (Vectorof Real) Real))
(define (vmin vec)
  (for/fold ([current : Real +inf.0])
            ([v (in-vector vec)])
    (if (< v current) v current)))

(: vmax (-> (Vectorof Real) Real))
(define (vmax vec)
  (for/fold ([current : Real -inf.0])
            ([v (in-vector vec)])
    (if (> v current) v current)))

(: vnormalize* (-> (Vectorof Real) (Vectorof Real)))
(define (vnormalize* vec)
  (define vec-min (vmin vec))
  (define vec-max (vmax vec))
  (v/ (vector-map (λ ([v : Real]) (- v vec-min)) vec) (- vec-max vec-min)))

(: vargmin (-> (Vectorof Real) Index))
(define (vargmin vec)
  (for/fold ([cur-val : Real +inf.0]
             [cur-idx : Integer 0]
             #:result (cast cur-idx Index))
            ([(v idx) (in-indexed (in-vector vec))])
    (if (< v cur-val)
        (values v idx)
        (values cur-val cur-idx))))

(: vargmax (-> (Vectorof Real) Index))
(define (vargmax vec)
  (for/fold ([cur-val : Real -inf.0]
             [cur-idx : Integer 0]
             #:result (cast cur-idx Index))
            ([(v idx) (in-indexed (in-vector vec))])
    (if (> v cur-val)
        (values v idx)
        (values cur-val cur-idx))))

(: tricubic (-> Real Nonnegative-Real))
(define (tricubic n)
  (cond [(> (abs n) 1) 0]
        [else (cast
               (expt (- 1 (expt (abs n) 3)) 3)
               Nonnegative-Real)]))

(: get-indexer (-> (Vectorof Real) Integer (Listof Integer)))
(define (get-indexer distances window)
  (define idx (vargmin distances))
  (define n (vector-length distances))

  (define left (max 0 (- idx (quotient window 2))))
  (define right
    (if (> (+ left window) (sub1 n))
        (sub1 n)
        (+ left window)))

  (range left (add1 right)))

(: index-with (All (A) (-> (Vectorof A) (Listof Index) (Vectorof A))))
(define (index-with vec indexer)
  (for/vector : (Vectorof A) ([idx (in-list indexer)])
    (vector-ref vec idx)))

(: get-weights (-> (Vectorof Real) (Listof Index) (Vectorof Real)))
(define (get-weights distances indexer)
  (define indexed (index-with distances indexer))
  (define max-dist (vmax indexed))
  (vector-map tricubic (v/ indexed max-dist)))

(: normalize-x (-> (Vectorof Real) Real Real))
(define (normalize-x xs x)
  (define xs-min (vmin xs))
  (define xs-max (vmax xs))
  (/ (- x xs-min) (- xs-max xs-min)))

(: denormalize-y (-> (Vectorof Real) Real Real))
(define (denormalize-y ys y)
  (define ys-min (vmin ys))
  (define ys-max (vmax ys))
  (+ (* y (- ys-max ys-min)) ys-min))

(: loess-fit (->* [(Vectorof Real) (Vectorof Real) Integer]
                  [#:degree Positive-Integer]
                  (-> Real Real)))
(define (loess-fit xs ys window #:degree [degree 1])
  (define xs-norm (vnormalize* xs))
  (define ys-norm (vnormalize* ys))
  (λ (x)
    (define norm-x (normalize-x xs x))
    (define distances (vector-map (λ ([v : Real]) (abs (- v norm-x))) xs-norm))
    (define indexer (cast (get-indexer distances window) (Listof Index)))
    (define weights (get-weights distances indexer))

    (define W (diagonal-matrix (vector->list weights)))
    ; XXX: gross
    (define indexed (index-with xs-norm indexer))
    (define X^T
      (vector*->matrix
       (for/vector : (Vectorof (Vectorof Real)) ([i (in-range (add1 degree))])
         (vector-map (λ ([v : Real]) (cast (expt v i) Real)) indexed))))
    (define X (matrix-transpose X^T))
    (define X^T*W (matrix* X^T W))
    (define Y (->col-matrix (index-with ys-norm indexer)))
    (define beta (matrix-transpose
                  (matrix* (matrix-inverse (matrix* X^T*W X))
                           X^T*W Y)))

    (define Xp (->col-matrix (build-list (add1 degree) (λ ([p : Number]) (cast (expt norm-x p) Real)))))
    (define normalized (matrix* beta Xp))
    (denormalize-y ys (matrix-ref normalized 0 0))))

(define xx : (Vectorof Real)
  (vector 0.5578196 2.0217271 2.5773252 3.4140288 4.3014084
          4.7448394 5.1073781 6.5411662 6.7216176 7.2600583
          8.1335874 9.1224379 11.9296663 12.3797674 13.2728619
          14.2767453 15.3731026 15.6476637 18.5605355 18.5866354
          18.7572812))
(define yy : (Vectorof Real)
  (vector 18.63654 103.49646 150.35391 190.51031 208.70115
          213.71135 228.49353 233.55387 234.55054 223.89225
          227.68339 223.91982 168.01999 164.95750 152.61107
          160.78742 168.55567 152.42658 221.70702 222.69040
          243.18828))

(plot (list (points (vector-map (inst vector Real) xx yy))
            (function (loess-fit xx yy 7 #:degree 1))))
