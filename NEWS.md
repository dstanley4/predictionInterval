# predictionInterval 1.0.1

* The MBESS package is no longer required. The confidence interval for d used by pi.d is now calculated within predictionInterval, using the same noncentral t approach as before. Prediction intervals are unchanged, with the exceptions listed below.

* The ggplot2 package is no longer required. It was imported but never used.

* pi.d Intervals are now correct when the noncentrality parameter is beyond 37.62 (a very large d-value, or a moderate d-value with thousands of participants per group in the original or replication study). In earlier versions a limit could be off in the second decimal in these cases, and by more for d-values far beyond the usual range.

* pi.d and pi.d.demo No longer turn off warnings for the rest of the R session, and no longer print a message about the noncentrality parameter exceeding 37.62.

* prob.level must now be a value with at most two decimals from .01 to .99 (for example .90 or .95), in every function. Any other value stops with a message. In earlier versions levels such as .999 were calculated but could not be printed, and pi.d returned a meaningless interval for a value such as 95.

* pi.d Now stops with a message when d is not a single finite number.

* pi.d Now works when cell sizes above 46,340 are given as integers.

* pi.d The interpretation text now refers to the prediction interval (PI). It previously called it a replication interval (RI).

* pi.d, pi.r and pi.m The interpretation text has been reworded (for example "based on a replication sample size"), and the printed output now ends with a line break.

* pi.d, pi.r and pi.m Printing no longer fails at some probability levels (for example .57).

* Sample sizes must be whole numbers. Other values stop with a message; previously they stopped with a message about an invalid format.

* pi.m.demo and pi.r.demo The illustrative trials now print the upper confidence limit in the ci.UL column. Previously the lower limit was printed in both columns (#2).

* pi.r.demo The note under the illustrative trials now labels the lower prediction limit pi.LL. It was labelled pi.UL.

* pi.d.demo The illustrative trials are now printed to two decimals, as in pi.r.demo and pi.m.demo.

* Corrections to the help pages: parameter descriptions that named the wrong statistic, typing errors, and the reference to Spence and Stanley (2016), which was listed as in preparation.

* pi.d.demo Runs about seven times faster in most cases. It is slower than before (up to about five times) when the noncentrality parameter is beyond 37.62, where the earlier version was fast but inaccurate.
