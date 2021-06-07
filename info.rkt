#lang info

(define version "0.1")
(define collection "loess")
(define name "loess")
(define deps '("base"
               "math-lib"
               "plot-lib"
               "typed-racket-lib"))

(define scribblings '(("scribblings/loess.scrbl" ())))
(define build-deps '("plot-doc"
                     "plot-gui-lib"
                     "racket-doc"
                     "scribble-lib"
                     "typed-racket-doc"))