#ifndef TRACK_SET_MACROS_HH
#define TRACK_SET_MACROS_HH

#define TRKCOV_FROMARRAY(PAR1,PAR2)			\
    cov(PAR1,PAR2) = cov_array[PAR1 ## PAR2];	\
    cov(PAR2,PAR1) = cov_array[PAR1 ## PAR2]
#define TRKCOV_2FROMARRAY(PAR) cov(PAR,PAR) = cov_array[PAR ## PAR]

#define TRKCOV_TOARRAY(PAR1,PAR2)			\
    cov_array[PAR1 ## PAR2] = cov(PAR1,PAR2);	\
    cov_array[PAR1 ## PAR2] = cov(PAR2,PAR1)
#define TRKCOV_2TOARRAY(PAR) cov_array[PAR ## PAR] = cov(PAR,PAR)

#endif
