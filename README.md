# br1
Contain sw 
These set of programs implement least square collocation method.
It pass througth several steps:

0) Detrend the date with polinomial or spline approximation 
1) Calculate optimun lag (Corona prog)
2) Calculate autocovariance and crosscoveriance empirical function (Correl)
3) Find best autocovariance crosscovariance model (Intcor)
4) Apply collocation filtering (collo-t/p/s)
5) Estrapolate componente taking into account solution of collocation
