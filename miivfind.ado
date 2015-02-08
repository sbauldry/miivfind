
capture program drop miivfind
*!miivfind v1.0.0 SBauldry 13sep2012

program define miivfind, rclass
	version 12.1
	syntax namelist
	
	*** constructing macros for matrix inputs
	tokenize `namelist'
	
	*** verifying basic vector/matrix entry
	* y1, y2, x1, x2 should be row vectors
	if (rowsof(`1') > 1) {
		dis as error `"Error: y1 is not a row vector."'
		exit
	}
	
	if (rowsof(`2') > 1) {
		dis as error `"Error: y2 is not a row vector."'
		exit
	}
	
	if (rowsof(`3') > 1) {
		dis as error `"Error: x1 is not a row vector."'
		exit
	}
	
	if (rowsof(`4') > 1) {
		dis as error `"Error: x2 is not a row vector."'
		exit
	}
	
	* ThetaE, ThetaD, and Psi should be square matrices
	if (rowsof(`9') != colsof(`9')) {
		dis as err `"Error: ThetaE is not a square matrix."'
		exit
	}
	
	if (rowsof(`10') != colsof(`10')) {
		dis as err `"Error: ThetaD is not a square matrix."'
		exit
	}
	
	if (rowsof(`11') != colsof(`11')) {
		dis as err `"Error: Psi is not a square matrix."'
		exit
	}
	
	*** verifying conformability of matrices
	* number of cols of y1 should equal number of cols of Beta
	if (colsof(`1') != colsof(`5')) {
		dis as err `"Error: y1 does not conform with Beta."'
		exit
	}
	
	* number of cols of y1 should equal number of rows of Gamma
	if (colsof(`1') != rowsof(`6')) {
		dis as err `"Error: y1 does not conform with Gamma."'
		exit
	}
	
	* number of cols of x1 should equal number of cols of Gamma
	if (colsof(`3') != colsof(`6')) {
		dis as err `"Error: x1 does not conform with Gamma."'
		exit
	}
	
	* number cols y2 should equal number rows Ly2
	if (colsof(`2') != rowsof(`7')) {
		dis as err `"Error: y2 does not conform with Ly2."'
		exit
	}
	
	* number cols x2 should equal number rows Lx2
	if (colsof(`4') != rowsof(`8')) {
		dis as err `"Error: x2 does not conform with Lx2."'
		exit
	}	
	
	
	*** running mata function to find MIIVs
	mata: miivfinder("`1'", "`2'", "`3'", "`4'", "`5'", "`6'", "`7'", "`8'", "`9'", "`10'", "`11'")
	
	*** construct header for output
	dis _newline(1)
	dis as text "List of MIIVs (if any) by Equation DV"
	dis as text "{hline 70}"
	dis as text "DV    MIIVs"
	dis as text "{hline 70}"

	*** displaying results in output
	local r = rowsof(IV)
	local c = colsof(IV)
	
	*** looping over each equation
	forval i = 1/`r' {
		
		*** finding location of first nonzero entry after column 1
		local sc = 0
		forval j  = 2/`c' {
			if IV[`i',`j'] != 0 {
				local sc = `j'
				continue, break
			}
		}

		*** counting number of MIIVs
		local nmiiv = 0
		forval j = 2/`c' {
			if IV[`i',`j'] != 0 {
				local nmiiv = `nmiiv' + 1
			}
		}

		*** preparing list of MIIVs
		if `nmiiv' == 0 local miiv = ""

		else if `nmiiv' == 1 local miiv = IV[`i',`sc']

		else if `nmiiv' > 1 & `sc' < `c' {
			local miiv = IV[`i',`sc']
			local ssc = `sc' + 1
			forval j = `ssc'/`c' {
				if IV[`i',`j'] != 0 {
					local miiv `miiv' ", " IV[`i',`j']
				}
			}
		}
		
		*** outputting row
		dis as text %3s IV[`i',1] "     " `miiv'
	}
	
	*** construct footer for output
	dis as text "{hline 70}"
	dis as text "Note: numbers in table are indices assigned to variables."
end


version 12.1
mata:
mata clear
void miivfinder(v1, v2, v3, v4, m1, m2, m3, m4, m5, m6, m7) 
{

// Defining vectors and matrices	
real rowvector y1
real rowvector y2
real rowvector x1
real rowvector x2
	
real matrix Beta
real matrix Gamma
real matrix Ly2
real matrix Lx2
	
real matrix ThetaE
real matrix ThetaD
real matrix Psi

// Reading vectors and matrices from Stata
y1 = st_matrix(v1)
y2 = st_matrix(v2)
x1 = st_matrix(v3)
x2 = st_matrix(v4)
	
Beta = st_matrix(m1)
Gamma = st_matrix(m2)
Ly2 = st_matrix(m3)
Lx2 = st_matrix(m4)
	
ThetaE = st_matrix(m5)
ThetaD = st_matrix(m6)
Psi = st_matrix(m7)	
	

// Generating predictors matrix (P)

// Starting with y1 equations	
if (y1 != (0)) {
	for (i = 1; i <= rows(Beta); i++) {
		Pi = y1[i]
		
		// Checking Beta matrix
		for (j = 1; j <= cols(Beta); j++) {
			if (Beta[i,j] !=0) Pi = Pi, y1[j]
		}
		
		// Checking Gamma matrix
		if (x1 != (0)) {
			for (j = 1; j <= cols(Gamma); j++) {
				if (Gamma[i,j] != 0) Pi = Pi, x1[j]
			}
		}
		
		// Ensuring conformability for stacking
		if (i == 1) P = Pi
		else {
			if ( cols(P) > cols(Pi) ) {
				add = J(1, ( cols(P) - cols(Pi) ), 0)
				Pi = Pi, add
			}
			else if ( cols(Pi) > cols(P) ) {
				add = J( rows(P), ( cols(Pi) - cols(P) ), 0 )
				P = P, add
			}
			P = P \ Pi
		}
	}
}

// Adding equations from y2	
if (y2 != (0)) {
	for (i = 1; i <= cols(y2); i++) {
		Pi = y2[i]

		// Checking Lambda_y2 matrix
		for (j = 1; j <= cols(y1); j++) {
			if (Ly2[i,j] != 0) Pi = Pi, y1[j]
		}
		
		// Ensuring conformability for stacking	
		if ( cols(P) == cols(Pi) ) P = P \ Pi
		else {
			if ( cols(P) > cols(Pi) ) {
				add = J(1, ( cols(P) - cols(Pi) ), 0)
				Pi = Pi, add
			}
			else if ( cols(Pi) > cols(P) ) {
				add = J( rows(P), ( cols(Pi) - cols(P) ), 0 )
				P = P, add
			}
			P = P \ Pi
		}
	}
}

// Adding equations from x2	
if (x2 != (0)) {
	for (i = 1; i <= cols(x2); i++) {
		Pi = x2[i]

		// Checking Lambda_x2 matrix
		for (j = 1; j <= cols(x1); j++) {
			if (Lx2[i,j] != 0) Pi = Pi, x1[j]
		}

		// Ensuring conformability for stacking	
		if ( cols(P) == cols(Pi) ) P = P \ Pi
		else {
			if ( cols(P) > cols(Pi) ) {
				add = J(1, ( cols(P) - cols(Pi) ), 0)
				Pi = Pi, add
			}
			else if ( cols(Pi) > cols(P) ) {
				add = J( rows(P), ( cols(Pi) - cols(P) ), 0 )
				P = P, add
			}
			P = P \ Pi
		}
	}
}



// Generating composite matrix (C)

// Stacking y1 and y2, sorting, and creating index for Theta
if (y1 != (0)) {
	Y = y1'
	if (y2 != (0)) Y = Y \ y2'
	Y = sort(Y,1)
	Yo = 1
	if (length(Y) > 1) {
		for (i = 2; i <= length(Y); i++) {
			Yo = Yo \ i
		}
	}
	Y = Y, Yo
}

// Stacking x1 and x2, sorting, and creating index for Theta
if (x1 != (0)) {
	X = x1'
	if (x2 != (0)) X = X \ x2'
	X = sort(X,1)
	Xo = 1
	if (length(X) > 1) {
		for (i = 2; i <= length(X); i++) {
			Xo = Xo \ i
		}
	}
	X = X, Xo
}

// Find composite disturbance for each equation
for (i = 1; i <= rows(P); i++) {
	Ci = P[i,1]
	
	// Find zetas for DVs in y1
	if (y1 != (0)) {
		for (j = 1; j <= cols(y1); j++) {
			if (Ci == y1[j]) Ci = Ci, Psi[j,j]
		}
	}

	// Find other disturbances in Y
	for (j = 1; j <= cols(P); j++) {
		if (Y != (0)) {
			for (k = 1; k <= rows(Y); k++) {
				if ( P[i,j] == Y[k,1] ) Ci = Ci, ThetaE[ Y[k,2], Y[k,2] ]
			}
		}
		
		if (X != (0)) {
			for (k = 1; k <= rows(X); k++) {
				if ( P[i,j] == X[k,1] ) Ci = Ci, ThetaD[ X[k,2], X[k,2] ]
			}
		}
	}
	
	// Ensuring conformability for stacking
	if (i == 1) C = Ci
	else if ( cols(C) == cols(Ci) ) C = C \ Ci
	else {
		if ( cols(C) > cols(Ci) ) {
			add = J(1, ( cols(C) - cols(Ci) ), 0)
			Ci = Ci, add
		}
		else if ( cols(C) < cols(Ci) ) {
			add = J( rows(C), ( cols(Ci) - cols(C) ), 0 )
			C = C, add
		}
		C = C \ Ci
	}

}	



// Generating total effects matrix (T)

// Total effects of errors on Y
if (y1 != (0)) {
	
	// total effects of epsilon on Y
	TotE = diagonal(ThetaE)
	
	// total effects of zeta on y1 through Beta
	TotZy1 = luinv( I( cols(y1) ) - Beta )
	
	// total effects of zeta on y2 through Ly2
	if (y2 != (0)) TotZy2 = Ly2*TotZy1
	
	// filling in index for zeta
	for (i = 1; i <= rows(TotZy1); i++) {
		for (j = 1; j <= cols(TotZy1); j++) {
			if (TotZy1[i,j] != 0) TotZy1[i,j] = Psi[j,j]
		}
	}
	
	if (y2 != (0)) {
		for (i = 1; i <= rows(TotZy2); i++) {
			for (j = 1; j <= cols(TotZy2); j++) {
				if (TotZy2[i,j] != 0) TotZy2[i,j] = Psi[j,j]
			}
		}
	}
	
	// add index, combine TotZy1 and TotZy2, and sort
	TotZy1 = y1', TotZy1
	if (y2 != (0)) {
		TotZy2 = y2', TotZy2
		TotY = TotZy1 \ TotZy2
	}
	else {
		TotY = TotZy1
	}
	
	TotY = TotE, sort(TotY,1)
	
	// setting T if no x1
	if (x1 == (0)) T = TotY
}

// Total effects of errors on X
if (x1 != (0)) {
	if (x2 != (0)) TotX = x1' \ x2'
	else TotX = x1'
	
	// total effects of delta on X
	TotD = diagonal(ThetaD)
	
	// sorting and adding
	TotX = sort(TotX,1), TotD
	
	// setting T if no y1
	if (y1 == (0)) T = TotX
}

// Combining TotY and TotX

// Ensuring conformable
if (y1 != (0) & x1 != (0)) {
	if ( cols(TotY) > cols(TotX) ) {
		add = J( rows(TotX), cols(TotY) - cols(TotX), 0 )
		TotX = TotX, add
	}
	else if ( cols(TotY) < cols(TotX) ) {
		add = J( rows(TotY), cols(TotX) - cols(TotY), 0 )
		TotY = TotY, add
	}
	T = TotY \ TotX
}


// Generating matrix of potential IVs (PIV)

for (i = 1; i <= rows(C); i++) {
	PIVi = C[i,1]
	for (p = 1; p <= rows(T); p++) {
		inst = 1	
		for (j = 2; j <= cols(C); j++) {
			if ( C[i,j] > 0 & inst == 1) {
				for (q = 2; q <= cols(T); q++) {
					if ( C[i,j] == T[p,q] ) inst = 0
				}
			}
		}
		if ( inst == 1 ) PIVi = PIVi, T[p,1]
	}
	
	// Ensuring conformability for stacking
	if (i == 1) PIV = PIVi
	else if ( cols(PIV) == cols(PIVi) ) PIV = PIV \ PIVi
	else {
		if ( cols(PIV) > cols(PIVi) ) {
			add = J(1, ( cols(PIV) - cols(PIVi) ), 0)
			PIVi = PIVi, add
		}
		else if ( cols(PIV) < cols(PIVi) ) {
			add = J( rows(PIV), ( cols(PIVi) - cols(PIV) ), 0 )
			PIV = PIV, add
		}
		PIV = PIV \ PIVi
	}
}


// Generating matrix of IVs (IV)

// Checking for correlated disturbances

IV = PIV

// Loop over equations and potential instruments in PIV
for (i = 1; i <= rows(PIV); i++) {
	for (p = 2; p <= cols(PIV); p++) {
	
		// if potential instrument
		if ( PIV[i,p] != 0 ) {
		
			// loop over disturbances in each equation
			for (j = 2; j <= cols(C); j++) {
			
				// loop over variables disturbances effect
				for (q = 2; q <= cols(T); q++) {
				
					// ???
					if ( C[i,j] != 0 & T[ PIV[i,p],q ] != 0 ) {
						
						// check for covariances in ThetaE
						if ( sum( ThetaE :== C[i,j] ) != 0 & sum( ThetaE :== T[ PIV[i,p],q ] ) != 0 ) {				
							for (r = 1; r <= rows(ThetaE); r++) {
								for (c = 1; c <= cols(ThetaE); c++) {
									if (C[i,j] == ThetaE[r,c]) loc1 = r
									if (T[ PIV[i,p],q ] == ThetaE[r,c]) loc2 = c
								}
							}
							if ( ThetaE[ loc1, loc2 ] != 0 ) IV[i,p] = 0
						}
						
						// check for covariances in ThetaD
						if ( sum( ThetaD :== C[i,j] ) != 0 & sum( ThetaD :== T[ PIV[i,p],q ] ) != 0 ) {				
							for (r = 1; r <= rows(ThetaD); r++) {
								for (c = 1; c <= cols(ThetaD); c++) {
									if (C[i,j] == ThetaD[r,c]) loc1 = r
									if (T[ PIV[i,p],q ] == ThetaD[r,c]) loc2 = c
								}
							}
							if ( ThetaD[ loc1, loc2 ] != 0 ) IV[i,p] = 0
						}
						
						// check for covariances in Psi
						if ( sum( Psi :== C[i,j] ) != 0 & sum( Psi :== T[ PIV[i,p],q ] ) != 0 ) {				
							for (r = 1; r <= rows(Psi); r++) {
								for (c = 1; c <= cols(Psi); c++) {
									if (C[i,j] == Psi[r,c]) loc1 = r
									if (T[ PIV[i,p],q ] == Psi[r,c]) loc2 = c
								}
							}
							if ( Psi[ loc1, loc2 ] != 0 ) IV[i,p] = 0
						}
					}
				}
			}
		}
	}
}

// Returning IV matrix to Stata
st_matrix("IV",IV)
}

end

