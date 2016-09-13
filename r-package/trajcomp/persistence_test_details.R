

alpha  <- c(0, -1, 5.5, 0, -1,2,4.8,3,3.5,4.2,4.3,9,3,0)


#alpha  <- c(0, 9,0 ,-7,-2,-3,6, -2, 1,0)

#minima <- persistance_minima(alpha)

debug_plot <- function(bars, comps, iter){
  bars_min <- bars[c(TRUE, FALSE)]
  bars_max <- bars[c(FALSE, TRUE)]
  
  comps_left <- comps[c(TRUE, FALSE, FALSE)]
  comps_min <- comps[c(FALSE, TRUE, FALSE)]
  comps_right <- comps[c(FALSE, FALSE, TRUE)]
  
  # curve and bar plot
  plot(alpha, pch=1, col="gray11", ylab="Curvature", xlab = "Index", main="Component growth")
  lines(alpha, col="gray11", lty=2)
  if(iter >= 0)
    text(12.5, -0.75, paste("iteration",as.character(iter)))
  #lines(alpha, col="orange")
  
  xLow <- bars_min+1
  yLow <- alpha[bars_min+1]
  
  xHigh <- bars_max+1
  yHigh <- alpha[bars_max+1]
  
  #arrows(xHigh,yHigh,xLow,yLow,col="orange",angle=90,length=0.1,code=3)
  
  
  comp_l_x <- comps_left+1
  comp_l_y <- alpha[comps_left+1]
  
  comp_r_x <- comps_right+1
  comp_r_y <- alpha[comps_right+1]
  
  comp_m_x <- comps_min+1
  comp_m_y <- alpha[comps_min+1]
  
  #arrows(comp_l_x,comp_l_y,comp_r_x,comp_r_y,col="purple", lty=3,angle=45,length=0.05,code=3, lwd=2)
  arrows(comp_l_x,comp_l_y,comp_m_x,comp_m_y,col="blue", lty=1,angle=45,length=0,code=3, lwd=2)
  arrows(comp_m_x,comp_m_y,comp_r_x,comp_r_y,col="red", lty=1,angle=45,length=0,code=3, lwd=2)
  
  points(bars_min+1, alpha[bars_min+1], col="green", pch=4, lwd=2)
  points(bars_max+1, alpha[bars_max+1], col="red", pch=2, lwd=2)
  legend('topleft',
         c("curve", "component left", "component right", "maxima", "minima"),
         col = c("gray11","blue", "red", "red", "green"),
         lty = c(2, 1, 1, -1,-1),
         pch = c(1, -1, -1, 2, 4),
         lwd=c(1,2,2,2,2),
         bty="n")
}

for(iter in 0:length(alpha)){
  bars <- persistence_test_bars(alpha, 0, -1)
  comps <- persistence_test_comps(alpha, 0, iter)
  debug_plot(bars,comps, iter)
  
  print(comps);
  print(bars);
}


bars <- persistence_test_bars(alpha, 0, -1)
comps <- persistence_test_comps(alpha, 0, -1)

debug_plot(bars,comps, -1)

