# TODO: Add comment
# 
# Author: zhaos
###############################################################################

uniroot.integer<-function (f, interval, ..., lower = min(interval), upper = max(interval), 
		step.power = 6, step.up = TRUE, pos.side = FALSE, print.steps = FALSE, 
		maxiter = 1000) 
{
	stored<-NULL
	iter <- 0
	if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
			upper) 
		stop("lower < upper  is not fulfilled")
	if (lower == -Inf && step.up == TRUE) 
		stop("lower cannot be -Inf when step.up=TRUE")
	if (upper == Inf && step.up == FALSE) 
		stop("upper cannot be Inf when step.up=FALSE")
	if (step.up) {
		f.old <- f(lower, ...)
		iter <- iter + 1
		sign <- 1
		xold <- lower
	}
	else {
		f.old <- f(upper, ...)
		iter <- iter + 1
		sign <- -1
		xold <- upper
	}
	if (print.steps) {
		print(paste("x=", xold, " f(x)=", f.old))
	}
	stored<-rbind(stored,c(xold,f.old))
	ever.switched <- FALSE
	tried.extreme <- FALSE
	while (step.power > -1) {
		xnew <- xold + sign * 2^step.power
		if ((step.up & xnew < upper) || (!step.up & xnew > lower)) {
			f.new <- f(xnew, ...)
			iter <- iter + 1
			if (print.steps) {
				print(paste("x=", xnew, " f(x)=", f.new))
			}
			stored<-rbind(stored,c(xnew,f.new))
		}
		else {
			xnew <- xold
			f.new <- f.old
			step.power <- step.power - 1
			if (tried.extreme == FALSE) {
				if (step.up) {
					f.extreme <- f(upper, ...)
					iter <- iter + 1
					x.extreme <- upper
				}
				else {
					f.extreme <- f(lower, ...)
					iter <- iter + 1
					x.extreme <- lower
				}
				tried.extreme <- TRUE
				xswitch <- x.extreme
				f.switch <- f.extreme
				if (print.steps) {
					print(paste("x=", x.extreme, " f(x)=", f.extreme))
				}
				stored<-rbind(stored,c(x.extreme,f.extreme))
				
				if (f.old>0 & f.extreme>0) {
					warning(paste0("Function result larger than 0 when lower and upper value used (",lower,": ",f.old,"; ",upper,": ",f.extreme,")\n"))
					return(list(root=upper))
				}
				if (f.old<0 & f.extreme<0) {
					warning(paste0("Function result less than 0 when lower and upper value used (",lower,": ",f.old,"; ",upper,": ",f.extreme,")\n"))
					return(list(root=lower))
				}
#				if (f.old * f.extreme >= 0) {
#					stop("f() at extremes not of opposite sign")
#				}
			}
		}
		if (f.old * f.new < 0) {
			sign <- sign * (-1)
			ever.switched <- TRUE
			xswitch <- xold
			f.switch <- f.old
		}
		if (ever.switched) {
			step.power <- step.power - 1
			if (step.power == -1) {
				(break)()
			}
		}
		xold <- xnew
		f.old <- f.new
	}
	if (pos.side) {
		root <- ifelse(f.new > 0, xnew, xswitch)
		f.root <- ifelse(f.new > 0, f.new, f.switch)
	}
	else {
		root <- ifelse(f.new < 0, xnew, xswitch)
		f.root <- ifelse(f.new < 0, f.new, f.switch)
	}
	colnames(stored)<-c("N","Power")
	return(list(iter = iter, f.root = f.root, root = root,process=stored))
}
