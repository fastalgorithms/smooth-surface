        subroutine numthetafour(numtets,nlams)
        integer  numtets(nlams),nlams 
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
         numtets( 1) =  3
         numtets( 2) =  1
      endif
      if(nlams.eq. 3) then
         numtets( 1) =  3
         numtets( 2) =  4
         numtets( 3) =  1
      endif
      if(nlams.eq. 4) then
         numtets( 1) =  3
         numtets( 2) =  7
         numtets( 3) =  5
         numtets( 4) =  1
      endif
      if(nlams.eq. 5) then
         numtets( 1) =  3
         numtets( 2) =  7
         numtets( 3) =  9
         numtets( 4) =  8
         numtets( 5) =  1
      endif
      if(nlams.eq. 6) then
         numtets( 1) =  3
         numtets( 2) =  7
         numtets( 3) = 10
         numtets( 4) = 12
         numtets( 5) =  8
         numtets( 6) =  1
      endif
      if(nlams.eq. 7) then
         numtets( 1) =  3
         numtets( 2) =  7
         numtets( 3) = 11
         numtets( 4) = 15
         numtets( 5) = 16
         numtets( 6) =  9
         numtets( 7) =  1
      endif
      if(nlams.eq. 8) then
         numtets( 1) =  3
         numtets( 2) =  7
         numtets( 3) = 11
         numtets( 4) = 15
         numtets( 5) = 19
         numtets( 6) = 19
         numtets( 7) =  7
         numtets( 8) =  1
      endif
      if(nlams.eq. 9) then
         numtets( 1) =  4
         numtets( 2) =  7
         numtets( 3) = 11
         numtets( 4) = 15
         numtets( 5) = 20
         numtets( 6) = 20
         numtets( 7) = 24
         numtets( 8) =  7
         numtets( 9) =  1
      endif
      if(nlams.eq.10) then
         numtets( 1) =  4
         numtets( 2) =  7
         numtets( 3) = 11
         numtets( 4) = 15
         numtets( 5) = 20
         numtets( 6) = 24
         numtets( 7) = 27
         numtets( 8) = 25
         numtets( 9) =  6
         numtets(10) =  1
      endif
      if(nlams.eq.11) then
         numtets( 1) =  4
         numtets( 2) =  7
         numtets( 3) = 11
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 24
         numtets( 7) = 28
         numtets( 8) = 32
         numtets( 9) = 29
         numtets(10) =  6
         numtets(11) =  1
      endif
      if(nlams.eq.12) then
         numtets( 1) =  4
         numtets( 2) =  8
         numtets( 3) = 11
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 24
         numtets( 7) = 28
         numtets( 8) = 33
         numtets( 9) = 35
         numtets(10) = 32
         numtets(11) =  5
         numtets(12) =  1
      endif
      if(nlams.eq.13) then
         numtets( 1) =  4
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 19
         numtets( 6) = 24
         numtets( 7) = 28
         numtets( 8) = 32
         numtets( 9) = 37
         numtets(10) = 38
         numtets(11) = 37
         numtets(12) =  6
         numtets(13) =  1
      endif
      if(nlams.eq.14) then
         numtets( 1) =  4
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 24
         numtets( 7) = 29
         numtets( 8) = 33
         numtets( 9) = 38
         numtets(10) = 40
         numtets(11) = 44
         numtets(12) = 35
         numtets(13) =  5
         numtets(14) =  1
      endif
      if(nlams.eq.15) then
         numtets( 1) =  4
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 33
         numtets( 9) = 38
         numtets(10) = 41
         numtets(11) = 45
         numtets(12) = 47
         numtets(13) = 40
         numtets(14) =  5
         numtets(15) =  1
      endif
      if(nlams.eq.16) then
         numtets( 1) =  5
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 33
         numtets( 9) = 37
         numtets(10) = 42
         numtets(11) = 47
         numtets(12) = 49
         numtets(13) = 49
         numtets(14) = 46
         numtets(15) =  6
         numtets(16) =  1
      endif
      if(nlams.eq.17) then
         numtets( 1) =  5
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 34
         numtets( 9) = 38
         numtets(10) = 43
         numtets(11) = 47
         numtets(12) = 51
         numtets(13) = 56
         numtets(14) = 56
         numtets(15) = 49
         numtets(16) =  5
         numtets(17) =  1
      endif
      if(nlams.eq.18) then
         numtets( 1) =  5
         numtets( 2) =  8
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 20
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 34
         numtets( 9) = 38
         numtets(10) = 43
         numtets(11) = 47
         numtets(12) = 51
         numtets(13) = 56
         numtets(14) = 59
         numtets(15) = 59
         numtets(16) = 51
         numtets(17) =  4
         numtets(18) =  1
      endif
      if(nlams.eq.19) then
         numtets( 1) =  5
         numtets( 2) =  9
         numtets( 3) = 12
         numtets( 4) = 16
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 34
         numtets( 9) = 38
         numtets(10) = 43
         numtets(11) = 47
         numtets(12) = 52
         numtets(13) = 56
         numtets(14) = 60
         numtets(15) = 65
         numtets(16) = 65
         numtets(17) = 53
         numtets(18) =  5
         numtets(19) =  1
      endif
      if(nlams.eq.20) then
         numtets( 1) =  5
         numtets( 2) =  9
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 29
         numtets( 8) = 34
         numtets( 9) = 38
         numtets(10) = 43
         numtets(11) = 47
         numtets(12) = 52
         numtets(13) = 56
         numtets(14) = 60
         numtets(15) = 64
         numtets(16) = 68
         numtets(17) = 66
         numtets(18) = 53
         numtets(19) =  4
         numtets(20) =  1
      endif
      if(nlams.eq.21) then
         numtets( 1) =  5
         numtets( 2) =  9
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 38
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 61
         numtets(15) = 65
         numtets(16) = 68
         numtets(17) = 71
         numtets(18) = 68
         numtets(19) = 54
         numtets(20) =  5
         numtets(21) =  1
      endif
      if(nlams.eq.22) then
         numtets( 1) =  6
         numtets( 2) =  9
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 61
         numtets(15) = 66
         numtets(16) = 68
         numtets(17) = 74
         numtets(18) = 74
         numtets(19) = 75
         numtets(20) = 63
         numtets(21) =  4
         numtets(22) =  1
      endif
      if(nlams.eq.23) then
         numtets( 1) =  6
         numtets( 2) =  9
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 61
         numtets(15) = 66
         numtets(16) = 70
         numtets(17) = 75
         numtets(18) = 78
         numtets(19) = 82
         numtets(20) = 78
         numtets(21) = 67
         numtets(22) =  4
         numtets(23) =  1
      endif
      if(nlams.eq.24) then
         numtets( 1) =  6
         numtets( 2) =  9
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 61
         numtets(15) = 66
         numtets(16) = 70
         numtets(17) = 75
         numtets(18) = 79
         numtets(19) = 82
         numtets(20) = 83
         numtets(21) = 83
         numtets(22) = 67
         numtets(23) =  5
         numtets(24) =  1
      endif
      if(nlams.eq.25) then
         numtets( 1) =  6
         numtets( 2) = 10
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 25
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 75
         numtets(18) = 79
         numtets(19) = 83
         numtets(20) = 88
         numtets(21) = 88
         numtets(22) = 89
         numtets(23) = 71
         numtets(24) =  4
         numtets(25) =  1
      endif
      if(nlams.eq.26) then
         numtets( 1) =  6
         numtets( 2) = 10
         numtets( 3) = 13
         numtets( 4) = 17
         numtets( 5) = 21
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 52
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 75
         numtets(18) = 80
         numtets(19) = 84
         numtets(20) = 88
         numtets(21) = 91
         numtets(22) = 92
         numtets(23) = 89
         numtets(24) = 69
         numtets(25) =  5
         numtets(26) =  1
      endif
      if(nlams.eq.27) then
         numtets( 1) =  6
         numtets( 2) = 10
         numtets( 3) = 14
         numtets( 4) = 17
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 75
         numtets(18) = 80
         numtets(19) = 84
         numtets(20) = 89
         numtets(21) = 92
         numtets(22) = 95
         numtets(23) = 97
         numtets(24) = 95
         numtets(25) = 75
         numtets(26) =  4
         numtets(27) =  1
      endif
      if(nlams.eq.28) then
         numtets( 1) =  6
         numtets( 2) = 10
         numtets( 3) = 14
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 34
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 80
         numtets(19) = 85
         numtets(20) = 89
         numtets(21) = 93
         numtets(22) = 98
         numtets(23) =102
         numtets(24) =102
         numtets(25) = 99
         numtets(26) = 80
         numtets(27) =  4
         numtets(28) =  1
      endif
      if(nlams.eq.29) then
         numtets( 1) =  6
         numtets( 2) = 10
         numtets( 3) = 14
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 80
         numtets(19) = 85
         numtets(20) = 89
         numtets(21) = 94
         numtets(22) = 98
         numtets(23) =102
         numtets(24) =103
         numtets(25) =108
         numtets(26) =101
         numtets(27) = 86
         numtets(28) =  4
         numtets(29) =  1
      endif
      if(nlams.eq.30) then
         numtets( 1) =  7
         numtets( 2) = 10
         numtets( 3) = 14
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 43
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 66
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 80
         numtets(19) = 85
         numtets(20) = 89
         numtets(21) = 94
         numtets(22) = 97
         numtets(23) =103
         numtets(24) =107
         numtets(25) =110
         numtets(26) =112
         numtets(27) =108
         numtets(28) = 84
         numtets(29) =  4
         numtets(30) =  1
      endif
      if(nlams.eq.31) then
         numtets( 1) =  7
         numtets( 2) = 10
         numtets( 3) = 14
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 80
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 94
         numtets(22) = 99
         numtets(23) =103
         numtets(24) =107
         numtets(25) =110
         numtets(26) =115
         numtets(27) =117
         numtets(28) =111
         numtets(29) = 85
         numtets(30) =  5
         numtets(31) =  1
      endif
      if(nlams.eq.32) then
         numtets( 1) =  7
         numtets( 2) = 11
         numtets( 3) = 14
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 30
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 80
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 94
         numtets(22) = 99
         numtets(23) =103
         numtets(24) =107
         numtets(25) =112
         numtets(26) =114
         numtets(27) =119
         numtets(28) =117
         numtets(29) =108
         numtets(30) = 89
         numtets(31) =  4
         numtets(32) =  1
      endif
      if(nlams.eq.33) then
         numtets( 1) =  7
         numtets( 2) = 11
         numtets( 3) = 15
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 94
         numtets(22) = 99
         numtets(23) =103
         numtets(24) =107
         numtets(25) =112
         numtets(26) =116
         numtets(27) =120
         numtets(28) =124
         numtets(29) =124
         numtets(30) =116
         numtets(31) = 93
         numtets(32) =  5
         numtets(33) =  1
      endif
      if(nlams.eq.34) then
         numtets( 1) =  7
         numtets( 2) = 11
         numtets( 3) = 15
         numtets( 4) = 18
         numtets( 5) = 22
         numtets( 6) = 26
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 94
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =108
         numtets(25) =112
         numtets(26) =117
         numtets(27) =121
         numtets(28) =125
         numtets(29) =128
         numtets(30) =131
         numtets(31) =127
         numtets(32) = 95
         numtets(33) =  4
         numtets(34) =  1
      endif
      if(nlams.eq.35) then
         numtets( 1) =  7
         numtets( 2) = 11
         numtets( 3) = 15
         numtets( 4) = 19
         numtets( 5) = 23
         numtets( 6) = 27
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 95
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =108
         numtets(25) =113
         numtets(26) =117
         numtets(27) =121
         numtets(28) =125
         numtets(29) =130
         numtets(30) =133
         numtets(31) =133
         numtets(32) =124
         numtets(33) =105
         numtets(34) =  4
         numtets(35) =  1
      endif
      if(nlams.eq.36) then
         numtets( 1) =  7
         numtets( 2) = 11
         numtets( 3) = 15
         numtets( 4) = 19
         numtets( 5) = 23
         numtets( 6) = 27
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 39
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 95
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =108
         numtets(25) =113
         numtets(26) =117
         numtets(27) =122
         numtets(28) =126
         numtets(29) =130
         numtets(30) =133
         numtets(31) =137
         numtets(32) =137
         numtets(33) =128
         numtets(34) =103
         numtets(35) =  4
         numtets(36) =  1
      endif
      if(nlams.eq.37) then
         numtets( 1) =  8
         numtets( 2) = 11
         numtets( 3) = 15
         numtets( 4) = 19
         numtets( 5) = 23
         numtets( 6) = 27
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 40
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 95
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =108
         numtets(25) =113
         numtets(26) =118
         numtets(27) =122
         numtets(28) =126
         numtets(29) =131
         numtets(30) =135
         numtets(31) =139
         numtets(32) =142
         numtets(33) =144
         numtets(34) =139
         numtets(35) =107
         numtets(36) =  4
         numtets(37) =  1
      endif
      if(nlams.eq.38) then
         numtets( 1) =  8
         numtets( 2) = 12
         numtets( 3) = 15
         numtets( 4) = 19
         numtets( 5) = 23
         numtets( 6) = 27
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 40
         numtets(10) = 44
         numtets(11) = 48
         numtets(12) = 53
         numtets(13) = 57
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 95
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =109
         numtets(25) =113
         numtets(26) =118
         numtets(27) =122
         numtets(28) =127
         numtets(29) =131
         numtets(30) =135
         numtets(31) =139
         numtets(32) =143
         numtets(33) =143
         numtets(34) =141
         numtets(35) =139
         numtets(36) =106
         numtets(37) =  4
         numtets(38) =  1
      endif
      if(nlams.eq.39) then
         numtets( 1) =  8
         numtets( 2) = 12
         numtets( 3) = 15
         numtets( 4) = 19
         numtets( 5) = 23
         numtets( 6) = 27
         numtets( 7) = 31
         numtets( 8) = 35
         numtets( 9) = 40
         numtets(10) = 44
         numtets(11) = 49
         numtets(12) = 53
         numtets(13) = 58
         numtets(14) = 62
         numtets(15) = 67
         numtets(16) = 71
         numtets(17) = 76
         numtets(18) = 81
         numtets(19) = 85
         numtets(20) = 90
         numtets(21) = 95
         numtets(22) = 99
         numtets(23) =104
         numtets(24) =109
         numtets(25) =113
         numtets(26) =118
         numtets(27) =122
         numtets(28) =127
         numtets(29) =131
         numtets(30) =136
         numtets(31) =140
         numtets(32) =142
         numtets(33) =147
         numtets(34) =151
         numtets(35) =151
         numtets(36) =139
         numtets(37) =111
         numtets(38) =  5
         numtets(39) =  1
      endif
      return
      end
