#lang typed/racket/base
(require math/matrix
         plot/utils
         racket/list
         racket/vector)
(provide loess-fit)

(: ∑ (-> (Vectorof Real) Real))
(define (∑ vec)
  (for/sum : Real ([v (in-vector vec)]) v))

(: vc* (-> (Vectorof Real) (Vectorof Real) (Vectorof Real)))
(define (vc* vec1 vec2)
  (for/vector : (Vectorof Real)
      ([v1 (in-vector vec1)]
       [v2 (in-vector vec2)])
    (* v1 v2)))

(: expt-real (-> Real Real Real))
(define (expt-real x b)
  (assert (expt x b) real?))

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
             #:result (assert cur-idx index?))
            ([(v idx) (in-indexed (in-vector vec))])
    (if (< v cur-val)
        (values v idx)
        (values cur-val cur-idx))))

(: vargmax (-> (Vectorof Real) Index))
(define (vargmax vec)
  (for/fold ([cur-val : Real -inf.0]
             [cur-idx : Integer 0]
             #:result (assert cur-idx index?))
            ([(v idx) (in-indexed (in-vector vec))])
    (if (> v cur-val)
        (values v idx)
        (values cur-val cur-idx))))

(: tricubic (-> Real Real))
(define (tricubic n)
  (cond [(> (abs n) 1) 0]
        [else (expt-real (- 1.0 (expt-real (abs n) 3.0)) 3.0)]))

(: get-indexer (-> (Vectorof Real) Integer (Listof Index)))
(define (get-indexer distances window)
  (: take* (∀ (A) (-> (Listof A) Integer (Listof A))))
  (define (take* lst n)
    (if (>= (length lst) n) (take lst n) lst))

  (define indexed (map (inst cons Real Index)
                       (vector->list distances)
                       (build-list (vector-length distances) (λ ([x : Index]) x))))
  (define sorted (sort indexed (λ ([a : (Pairof Real Index)]
                                   [b : (Pairof Real Index)])
                                 (< (car a) (car b)))))
  (define nearest (map (inst cdr Real Index) (take* sorted window)))
  (range (apply min nearest) (assert (add1 (apply max nearest)) index?)))

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

(: loess-fit (->* [(Vectorof Real) (Vectorof Real)]
                  [#:span Positive-Real
                   #:degree Positive-Integer]
                  (-> Real Real)))
(define (loess-fit xs ys #:span [span 0.75] #:degree [degree 1])
  (define n (vector-length xs))
  (when (not (= n (vector-length ys)))
    (error 'loess-fit "input vectors are not the same length: xs ~a, ys ~a"
           (vector-length xs) (vector-length ys)))
  (when (not (<= (/ (add1 degree) n) 1))
    (error 'loess-fit "span is not between (degree + 1)/n and 1: degree ~a, span ~a, n ~a"
           degree span n))

  (define window (assert (inexact->exact (ceiling (* span n))) integer?))
  (define xs-norm (vnormalize* xs))
  (define ys-norm (vnormalize* ys))

  (λ (x)
    (define norm-x (normalize-x xs x))
    (define distances (vector-map (λ ([v : Real]) (abs (- v norm-x))) xs-norm))
    (define indexer (get-indexer distances window))
    (define weights (get-weights distances indexer))

    (define (with-matrix)
      (define W (diagonal-matrix (vector->list weights)))
      ; XXX: gross
      (define indexed (index-with xs-norm indexer))
      (define X^T
        (vector*->matrix
         (for/vector : (Vectorof (Vectorof Real)) ([i (in-range (add1 degree))])
           (vector-map (λ ([v : Real]) (expt-real v i)) indexed))))
      (define X (matrix-transpose X^T))
      (define X^T*W (matrix* X^T W))
      (define Y (->col-matrix (index-with ys-norm indexer)))

      (define β (matrix->vector
                 (matrix* (matrix-inverse (matrix* X^T*W X))
                          X^T*W Y)))
      (define Xp (build-vector (add1 degree) (λ ([p : Real]) (expt-real norm-x p))))
      (vdot β Xp))

    (define (with-linear)
      ; textbook weighted linear regression -- possibly faster
      (define indexed-x (index-with xs-norm indexer))
      (define indexed-y (index-with ys-norm indexer))

      (define ∑w (∑ weights))
      (define ∑wx (vdot indexed-x weights))
      (define ∑wy (vdot indexed-y weights))
      (define ∑wx^2 (vdot (vc* indexed-x indexed-x) weights))
      (define ∑wxy (vdot (vc* indexed-x indexed-y) weights))

      (define μx (/ ∑wx ∑w))
      (define μy (/ ∑wy ∑w))

      (define b (/ (- ∑wxy (* μx μy ∑w))
                   (- ∑wx^2 (* μx μx ∑w))))
      (define a (- μy (* b μx)))

      (+ a (* b norm-x)))

    (define normalized (if (= degree 1) (with-linear) (with-matrix)))
    (denormalize-y ys normalized)))
