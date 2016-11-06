% 	program pressure

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,4);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\Shu\exact_solution\Shu_p_exact.dat
    xx=Shu_p_exact(:,1);
    yy=Shu_p_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-5.0 5.0 0.7 12.3]);
   