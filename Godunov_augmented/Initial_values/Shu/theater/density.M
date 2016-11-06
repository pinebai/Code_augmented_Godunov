% 	program density

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,2);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\Shu\exact_solution\Shu_d_exact.dat
    xx=Shu_d_exact(:,1);
    yy=Shu_d_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-5.0 5.0 0.7 5.21]);
   