# LOWESS-in-R
LOWESS smoother in R

I've seen people use this smoother a lot and was curious how it worked and programmed it in R. It is basically a rolling regression.

![chart of explanation of algo](https://github.com/frogger21/LOWESS-in-R/blob/master/lowess1.png)

A more detailed explanation on: https://www.itl.nist.gov/div898/handbook/pmd/section1/pmd144.htm

Here is an example from NIST using a rolling local window of 7 observations.
![Example of LOWESS from NIST](https://github.com/frogger21/LOWESS-in-R/blob/master/lowess2.png)

Using more iterations and a higher order polynomial it looks a bit smoother.
![Example of LOWESS from NIST more smooth](https://github.com/frogger21/LOWESS-in-R/blob/master/lowess3.png)

Now an example of applying the LOWESS smoother to time series data: sunspot data from R’s internal database.  Using a 12 month window we capture the faster cycles that occur frequently while using a longer window of 480 months (40 years) we can gauge the longer term trends of sunspot activity. It’s important to note that the smoother is impacted by future data points so what the smoother shows near the end is something I wouldn’t pay much attention to. 
![R Sunspots Monthly](https://github.com/frogger21/LOWESS-in-R/blob/master/lowess4.png)

Here’s an example emphasizing why the output of the smoother near the ends can be misleading. The dashed line is the output of the LOWESS smoother estimated with only 4/5th of the entire time series sample. 
![R Sunspots Monthly misleading ends](https://github.com/frogger21/LOWESS-in-R/blob/master/lowess5.png)


7/24/2020 update
I added a density parameter to estimate the LOWESS value of unobserved variables. For instance on the example given by NIST https://www.itl.nist.gov/div898/handbook/pmd/section1/dep/dep144.htm they show that the unobserved x=10 has the LOWESS value of f(x)=202.9876. It is calculated in the same straightforward way as the observed variables. We look for the n-window closest points to the point x=10 and then run through the distance, scaled distance, tricube weights and then use the WLS regression to find the LOWESS value of x=10. 
