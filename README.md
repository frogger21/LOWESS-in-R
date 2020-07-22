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
