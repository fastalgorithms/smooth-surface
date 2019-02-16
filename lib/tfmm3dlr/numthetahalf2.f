	subroutine numthetahalf(numtets,nlams)
	integer *4 numtets(nlams),nlams 
c
c This routine returns the number of Fourier modes needed in the
c phi integral for each of the discrete lambda values given
c by Norman's quadratures. (see weights.f).
c
c Input arguments:
c	nlams - number of nodes in the lambda quadrature.  This must
c	  be an integer in the range [2,39].
c
c Output arguments:
c	numtest(i) - number of Fourier modes needed for phi
c                    integral with lambda_i.
c
c Approximate accuracies of the quadrature:
c
c	NLAMS	maximum error
c
c        2     0.15318E+00
c        3     0.76505E-01
c        4     0.32149E-01
c        5     0.15630E-01
c        6     0.75110E-02
c        7     0.35030E-02
c        8     0.16243E-02
c        9     0.72230E-03
c       10     0.33074E-03
c       11     0.15035E-03
c       12     0.70952E-04
c       13     0.31751E-04
c       14     0.14589E-04
c       15     0.64300E-05
c       16     0.29477E-05
c       17     0.13222E-05
c       18     0.61488E-06
c       19     0.27435E-06
c       20     0.12534E-06
c       21     0.55324E-07
c       22     0.25257E-07
c       23     0.11293E-07
c       24     0.52063E-08
c       25     0.23256E-08
c       26     0.10580E-08
c       27     0.46835E-09
c       28     0.21286E-09
c       29     0.95164E-10
c       30     0.43599E-10
c       31     0.19516E-10
c       32     0.88491E-11
c       33     0.39313E-11
c       34     0.17821E-11
c       35     0.79603E-12
c       36     0.36460E-12
c       37     0.16331E-12
c       38     0.73497E-13
c       39     0.31530E-13
c
      if(nlams.eq. 2) then
         numtets( 1) = 1
         numtets( 2) = 1
      endif
      if(nlams.eq. 3) then
         numtets( 1) = 1
         numtets( 2) = 2
         numtets( 3) = 2
      endif
      if(nlams.eq. 4) then
         numtets( 1) = 1
         numtets( 2) = 2
         numtets( 3) = 2
         numtets( 4) = 1
      endif
      if(nlams.eq. 5) then
         numtets( 1) = 2
         numtets( 2) = 2
         numtets( 3) = 3
         numtets( 4) = 3
         numtets( 5) = 1
      endif
      if(nlams.eq. 6) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 3
         numtets( 4) = 3
         numtets( 5) = 2
         numtets( 6) = 1
      endif
      if(nlams.eq. 7) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 3
         numtets( 4) = 4
         numtets( 5) = 4
         numtets( 6) = 3
         numtets( 7) = 1
      endif
      if(nlams.eq. 8) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 4
         numtets( 4) = 4
         numtets( 5) = 5
         numtets( 6) = 5
         numtets( 7) = 5
         numtets( 8) = 1
      endif
      if(nlams.eq. 9) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 4
         numtets( 4) = 5
         numtets( 5) = 5
         numtets( 6) = 4
         numtets( 7) = 5
         numtets( 8) = 4
         numtets( 9) = 1
      endif
      if(nlams.eq.10) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 4
         numtets( 4) = 5
         numtets( 5) = 5
         numtets( 6) = 6
         numtets( 7) = 6
         numtets( 8) = 6
         numtets( 9) = 6
         numtets(10) = 1
      endif
      if(nlams.eq.11) then
         numtets( 1) = 2
         numtets( 2) = 3
         numtets( 3) = 4
         numtets( 4) = 5
         numtets( 5) = 6
         numtets( 6) = 6
         numtets( 7) = 6
         numtets( 8) = 7
         numtets( 9) = 7
         numtets(10) = 5
         numtets(11) = 1
      endif
      if(nlams.eq.12) then
         numtets( 1) = 2
         numtets( 2) = 4
         numtets( 3) = 5
         numtets( 4) = 5
         numtets( 5) = 6
         numtets( 6) = 6
         numtets( 7) = 7
         numtets( 8) = 7
         numtets( 9) = 7
         numtets(10) = 7
         numtets(11) = 5
         numtets(12) = 1
      endif
      if(nlams.eq.13) then
         numtets( 1) = 3
         numtets( 2) = 4
         numtets( 3) = 5
         numtets( 4) = 6
         numtets( 5) = 6
         numtets( 6) = 7
         numtets( 7) = 7
         numtets( 8) = 8
         numtets( 9) = 8
         numtets(10) = 8
         numtets(11) = 7
         numtets(12) = 5
         numtets(13) = 1
      endif
      if(nlams.eq.14) then
         numtets( 1) = 3
         numtets( 2) = 4
         numtets( 3) = 5
         numtets( 4) = 6
         numtets( 5) = 7
         numtets( 6) = 7
         numtets( 7) = 8
         numtets( 8) = 8
         numtets( 9) = 8
         numtets(10) = 9
         numtets(11) = 8
         numtets(12) = 8
         numtets(13) = 7
         numtets(14) = 1
      endif
      if(nlams.eq.15) then
         numtets( 1) = 3
         numtets( 2) = 4
         numtets( 3) = 5
         numtets( 4) = 6
         numtets( 5) = 7
         numtets( 6) = 7
         numtets( 7) = 8
         numtets( 8) = 8
         numtets( 9) = 9
         numtets(10) = 9
         numtets(11) = 9
         numtets(12) = 9
         numtets(13) = 8
         numtets(14) = 7
         numtets(15) = 1
      endif
      if(nlams.eq.16) then
         numtets( 1) = 3
         numtets( 2) = 4
         numtets( 3) = 5
         numtets( 4) = 6
         numtets( 5) = 7
         numtets( 6) = 8
         numtets( 7) = 8
         numtets( 8) = 9
         numtets( 9) = 9
         numtets(10) =10
         numtets(11) =10
         numtets(12) =10
         numtets(13) =10
         numtets(14) = 9
         numtets(15) = 8
         numtets(16) = 1
      endif
      if(nlams.eq.17) then
         numtets( 1) = 3
         numtets( 2) = 4
         numtets( 3) = 6
         numtets( 4) = 6
         numtets( 5) = 7
         numtets( 6) = 8
         numtets( 7) = 9
         numtets( 8) = 9
         numtets( 9) =10
         numtets(10) =10
         numtets(11) =10
         numtets(12) =10
         numtets(13) =11
         numtets(14) =10
         numtets(15) = 9
         numtets(16) = 8
         numtets(17) = 1
      endif
      if(nlams.eq.18) then
         numtets( 1) = 3
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 7
         numtets( 5) = 7
         numtets( 6) = 8
         numtets( 7) = 9
         numtets( 8) = 9
         numtets( 9) =10
         numtets(10) =10
         numtets(11) =11
         numtets(12) =11
         numtets(13) =11
         numtets(14) =11
         numtets(15) =11
         numtets(16) =11
         numtets(17) = 7
         numtets(18) = 1
      endif
      if(nlams.eq.19) then
         numtets( 1) = 3
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 7
         numtets( 5) = 8
         numtets( 6) = 8
         numtets( 7) = 9
         numtets( 8) =10
         numtets( 9) =10
         numtets(10) =11
         numtets(11) =11
         numtets(12) =12
         numtets(13) =12
         numtets(14) =12
         numtets(15) =12
         numtets(16) =12
         numtets(17) = 9
         numtets(18) = 9
         numtets(19) = 1
      endif
      if(nlams.eq.20) then
         numtets( 1) = 3
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 7
         numtets( 5) = 8
         numtets( 6) = 9
         numtets( 7) = 9
         numtets( 8) =10
         numtets( 9) =11
         numtets(10) =11
         numtets(11) =12
         numtets(12) =12
         numtets(13) =12
         numtets(14) =12
         numtets(15) =13
         numtets(16) =13
         numtets(17) =12
         numtets(18) =12
         numtets(19) = 8
         numtets(20) = 1
      endif
      if(nlams.eq.21) then
         numtets( 1) = 4
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 7
         numtets( 5) = 8
         numtets( 6) = 9
         numtets( 7) =10
         numtets( 8) =10
         numtets( 9) =11
         numtets(10) =12
         numtets(11) =12
         numtets(12) =12
         numtets(13) =13
         numtets(14) =13
         numtets(15) =13
         numtets(16) =13
         numtets(17) =13
         numtets(18) =13
         numtets(19) =12
         numtets(20) =10
         numtets(21) = 1
      endif
      if(nlams.eq.22) then
         numtets( 1) = 4
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 7
         numtets( 5) = 8
         numtets( 6) = 9
         numtets( 7) =10
         numtets( 8) =11
         numtets( 9) =11
         numtets(10) =12
         numtets(11) =12
         numtets(12) =13
         numtets(13) =13
         numtets(14) =13
         numtets(15) =14
         numtets(16) =14
         numtets(17) =14
         numtets(18) =14
         numtets(19) =14
         numtets(20) =12
         numtets(21) =10
         numtets(22) = 1
      endif
      if(nlams.eq.23) then
         numtets( 1) = 4
         numtets( 2) = 5
         numtets( 3) = 6
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) = 9
         numtets( 7) =10
         numtets( 8) =11
         numtets( 9) =12
         numtets(10) =12
         numtets(11) =13
         numtets(12) =13
         numtets(13) =14
         numtets(14) =14
         numtets(15) =14
         numtets(16) =14
         numtets(17) =15
         numtets(18) =15
         numtets(19) =14
         numtets(20) =14
         numtets(21) =14
         numtets(22) = 9
         numtets(23) = 1
      endif
      if(nlams.eq.24) then
         numtets( 1) = 4
         numtets( 2) = 5
         numtets( 3) = 7
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) =10
         numtets( 7) =10
         numtets( 8) =11
         numtets( 9) =12
         numtets(10) =13
         numtets(11) =13
         numtets(12) =14
         numtets(13) =14
         numtets(14) =14
         numtets(15) =15
         numtets(16) =15
         numtets(17) =15
         numtets(18) =15
         numtets(19) =15
         numtets(20) =15
         numtets(21) =15
         numtets(22) =14
         numtets(23) =10
         numtets(24) = 1
      endif
      if(nlams.eq.25) then
         numtets( 1) = 4
         numtets( 2) = 6
         numtets( 3) = 7
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) =10
         numtets( 7) =11
         numtets( 8) =11
         numtets( 9) =12
         numtets(10) =13
         numtets(11) =13
         numtets(12) =14
         numtets(13) =14
         numtets(14) =15
         numtets(15) =15
         numtets(16) =15
         numtets(17) =16
         numtets(18) =16
         numtets(19) =16
         numtets(20) =16
         numtets(21) =16
         numtets(22) =13
         numtets(23) =14
         numtets(24) =11
         numtets(25) = 1
      endif
      if(nlams.eq.26) then
         numtets( 1) = 4
         numtets( 2) = 6
         numtets( 3) = 7
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) =10
         numtets( 7) =11
         numtets( 8) =12
         numtets( 9) =12
         numtets(10) =13
         numtets(11) =14
         numtets(12) =14
         numtets(13) =15
         numtets(14) =15
         numtets(15) =16
         numtets(16) =16
         numtets(17) =16
         numtets(18) =16
         numtets(19) =17
         numtets(20) =17
         numtets(21) =17
         numtets(22) =16
         numtets(23) =15
         numtets(24) =14
         numtets(25) =11
         numtets(26) = 1
      endif
      if(nlams.eq.27) then
         numtets( 1) = 4
         numtets( 2) = 6
         numtets( 3) = 7
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) =10
         numtets( 7) =11
         numtets( 8) =12
         numtets( 9) =13
         numtets(10) =13
         numtets(11) =14
         numtets(12) =15
         numtets(13) =15
         numtets(14) =16
         numtets(15) =16
         numtets(16) =16
         numtets(17) =17
         numtets(18) =17
         numtets(19) =17
         numtets(20) =17
         numtets(21) =17
         numtets(22) =17
         numtets(23) =17
         numtets(24) =16
         numtets(25) =15
         numtets(26) =12
         numtets(27) = 1
      endif
      if(nlams.eq.28) then
         numtets( 1) = 4
         numtets( 2) = 6
         numtets( 3) = 7
         numtets( 4) = 8
         numtets( 5) = 9
         numtets( 6) =10
         numtets( 7) =11
         numtets( 8) =12
         numtets( 9) =13
         numtets(10) =14
         numtets(11) =14
         numtets(12) =15
         numtets(13) =15
         numtets(14) =16
         numtets(15) =16
         numtets(16) =17
         numtets(17) =17
         numtets(18) =17
         numtets(19) =18
         numtets(20) =18
         numtets(21) =18
         numtets(22) =18
         numtets(23) =18
         numtets(24) =18
         numtets(25) =17
         numtets(26) =15
         numtets(27) =12
         numtets(28) = 1
      endif
      if(nlams.eq.29) then
         numtets( 1) = 4
         numtets( 2) = 6
         numtets( 3) = 7
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =11
         numtets( 7) =12
         numtets( 8) =12
         numtets( 9) =13
         numtets(10) =14
         numtets(11) =15
         numtets(12) =15
         numtets(13) =16
         numtets(14) =16
         numtets(15) =17
         numtets(16) =17
         numtets(17) =18
         numtets(18) =18
         numtets(19) =18
         numtets(20) =18
         numtets(21) =19
         numtets(22) =19
         numtets(23) =19
         numtets(24) =19
         numtets(25) =18
         numtets(26) =16
         numtets(27) =17
         numtets(28) =11
         numtets(29) = 1
      endif
      if(nlams.eq.30) then
         numtets( 1) = 5
         numtets( 2) = 6
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =11
         numtets( 7) =12
         numtets( 8) =13
         numtets( 9) =13
         numtets(10) =14
         numtets(11) =15
         numtets(12) =16
         numtets(13) =16
         numtets(14) =17
         numtets(15) =17
         numtets(16) =18
         numtets(17) =18
         numtets(18) =18
         numtets(19) =19
         numtets(20) =19
         numtets(21) =19
         numtets(22) =19
         numtets(23) =19
         numtets(24) =19
         numtets(25) =19
         numtets(26) =19
         numtets(27) =19
         numtets(28) =15
         numtets(29) =12
         numtets(30) = 1
      endif
      if(nlams.eq.31) then
         numtets( 1) = 5
         numtets( 2) = 6
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =11
         numtets( 7) =12
         numtets( 8) =13
         numtets( 9) =14
         numtets(10) =14
         numtets(11) =15
         numtets(12) =16
         numtets(13) =16
         numtets(14) =17
         numtets(15) =18
         numtets(16) =18
         numtets(17) =18
         numtets(18) =19
         numtets(19) =19
         numtets(20) =19
         numtets(21) =20
         numtets(22) =20
         numtets(23) =20
         numtets(24) =20
         numtets(25) =20
         numtets(26) =20
         numtets(27) =20
         numtets(28) =18
         numtets(29) =19
         numtets(30) =11
         numtets(31) = 1
      endif
      if(nlams.eq.32) then
         numtets( 1) = 5
         numtets( 2) = 6
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =11
         numtets( 7) =12
         numtets( 8) =13
         numtets( 9) =14
         numtets(10) =15
         numtets(11) =15
         numtets(12) =16
         numtets(13) =17
         numtets(14) =17
         numtets(15) =18
         numtets(16) =18
         numtets(17) =19
         numtets(18) =19
         numtets(19) =20
         numtets(20) =20
         numtets(21) =20
         numtets(22) =20
         numtets(23) =21
         numtets(24) =21
         numtets(25) =21
         numtets(26) =21
         numtets(27) =20
         numtets(28) =20
         numtets(29) =19
         numtets(30) =17
         numtets(31) =13
         numtets(32) = 1
      endif
      if(nlams.eq.33) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =11
         numtets( 7) =12
         numtets( 8) =13
         numtets( 9) =14
         numtets(10) =15
         numtets(11) =16
         numtets(12) =16
         numtets(13) =17
         numtets(14) =18
         numtets(15) =18
         numtets(16) =19
         numtets(17) =19
         numtets(18) =20
         numtets(19) =20
         numtets(20) =20
         numtets(21) =21
         numtets(22) =21
         numtets(23) =21
         numtets(24) =21
         numtets(25) =21
         numtets(26) =21
         numtets(27) =20
         numtets(28) =21
         numtets(29) =19
         numtets(30) =19
         numtets(31) =18
         numtets(32) =14
         numtets(33) = 1
      endif
      if(nlams.eq.34) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =10
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =13
         numtets( 9) =14
         numtets(10) =15
         numtets(11) =16
         numtets(12) =17
         numtets(13) =17
         numtets(14) =18
         numtets(15) =19
         numtets(16) =19
         numtets(17) =20
         numtets(18) =20
         numtets(19) =20
         numtets(20) =21
         numtets(21) =21
         numtets(22) =21
         numtets(23) =22
         numtets(24) =22
         numtets(25) =22
         numtets(26) =22
         numtets(27) =22
         numtets(28) =22
         numtets(29) =22
         numtets(30) =21
         numtets(31) =21
         numtets(32) =18
         numtets(33) =13
         numtets(34) = 1
      endif
      if(nlams.eq.35) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 8
         numtets( 4) = 9
         numtets( 5) =11
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =14
         numtets( 9) =15
         numtets(10) =15
         numtets(11) =16
         numtets(12) =17
         numtets(13) =18
         numtets(14) =18
         numtets(15) =19
         numtets(16) =19
         numtets(17) =20
         numtets(18) =20
         numtets(19) =21
         numtets(20) =21
         numtets(21) =22
         numtets(22) =22
         numtets(23) =22
         numtets(24) =22
         numtets(25) =23
         numtets(26) =23
         numtets(27) =23
         numtets(28) =23
         numtets(29) =23
         numtets(30) =22
         numtets(31) =22
         numtets(32) =21
         numtets(33) =19
         numtets(34) =13
         numtets(35) = 1
      endif
      if(nlams.eq.36) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 8
         numtets( 4) =10
         numtets( 5) =11
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =14
         numtets( 9) =15
         numtets(10) =16
         numtets(11) =16
         numtets(12) =17
         numtets(13) =18
         numtets(14) =19
         numtets(15) =19
         numtets(16) =20
         numtets(17) =20
         numtets(18) =21
         numtets(19) =21
         numtets(20) =22
         numtets(21) =22
         numtets(22) =22
         numtets(23) =23
         numtets(24) =23
         numtets(25) =23
         numtets(26) =23
         numtets(27) =23
         numtets(28) =23
         numtets(29) =23
         numtets(30) =23
         numtets(31) =23
         numtets(32) =22
         numtets(33) =21
         numtets(34) =18
         numtets(35) =15
         numtets(36) = 1
      endif
      if(nlams.eq.37) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 8
         numtets( 4) =10
         numtets( 5) =11
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =14
         numtets( 9) =15
         numtets(10) =16
         numtets(11) =17
         numtets(12) =17
         numtets(13) =18
         numtets(14) =19
         numtets(15) =19
         numtets(16) =20
         numtets(17) =21
         numtets(18) =21
         numtets(19) =22
         numtets(20) =22
         numtets(21) =22
         numtets(22) =23
         numtets(23) =23
         numtets(24) =23
         numtets(25) =24
         numtets(26) =24
         numtets(27) =24
         numtets(28) =24
         numtets(29) =24
         numtets(30) =24
         numtets(31) =24
         numtets(32) =24
         numtets(33) =22
         numtets(34) =22
         numtets(35) =21
         numtets(36) =14
         numtets(37) = 1
      endif
      if(nlams.eq.38) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 9
         numtets( 4) =10
         numtets( 5) =11
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =14
         numtets( 9) =15
         numtets(10) =16
         numtets(11) =17
         numtets(12) =18
         numtets(13) =18
         numtets(14) =19
         numtets(15) =20
         numtets(16) =20
         numtets(17) =21
         numtets(18) =22
         numtets(19) =22
         numtets(20) =22
         numtets(21) =23
         numtets(22) =23
         numtets(23) =24
         numtets(24) =24
         numtets(25) =24
         numtets(26) =24
         numtets(27) =25
         numtets(28) =25
         numtets(29) =25
         numtets(30) =25
         numtets(31) =25
         numtets(32) =25
         numtets(33) =24
         numtets(34) =24
         numtets(35) =22
         numtets(36) =21
         numtets(37) =15
         numtets(38) = 1
      endif
      if(nlams.eq.39) then
         numtets( 1) = 5
         numtets( 2) = 7
         numtets( 3) = 9
         numtets( 4) =10
         numtets( 5) =11
         numtets( 6) =12
         numtets( 7) =13
         numtets( 8) =14
         numtets( 9) =15
         numtets(10) =16
         numtets(11) =17
         numtets(12) =18
         numtets(13) =19
         numtets(14) =19
         numtets(15) =20
         numtets(16) =21
         numtets(17) =21
         numtets(18) =22
         numtets(19) =22
         numtets(20) =23
         numtets(21) =23
         numtets(22) =24
         numtets(23) =24
         numtets(24) =24
         numtets(25) =25
         numtets(26) =25
         numtets(27) =25
         numtets(28) =25
         numtets(29) =25
         numtets(30) =26
         numtets(31) =25
         numtets(32) =25
         numtets(33) =25
         numtets(34) =25
         numtets(35) =25
         numtets(36) =23
         numtets(37) =20
         numtets(38) =15
         numtets(39) = 1
      endif
      return
      end
