// This do file performs a cubic spline interpolation in Stata, using a sample data// file posted on the-idea-shop.com. The file has 36 quarterly figures for Personal // income from the U.S. Bureau of Economics Analysis, which we want to interpolate// into monthly figures. //// See http://chamberlaineconomics.com/2010/01/20/how-economists-convert-quarterly-data-into-monthly-cubic-spline-interpolation///// (c) 2011 Chamberlain Economics, L.L.C. (info@chamberlaineconomics.com).use "https://www.dropbox.com/s/3y3h39dwbgghusp/stata_spline.dta?dl=0"mata // This line launches the mata system inside StataX = st_data((1,11),"x") // This pulls in the x quarterly markers data.Y = st_data((1,11),"y") // This pulls in the quarterly y data we want to interpolate between.XX = st_data(.,"xx") // This pulls in the xx monthly markers we want to interpolate at.A = spline3(X,Y) // This generates the cubic spline coefficients matrix, and stores it in A. B = spline3eval(A,XX) // This performs the interpolation, and store the values in B.st_store(.,"yy",B) // This pushes the inpolated figures in B back into the yy variable in Stata.endbrowse