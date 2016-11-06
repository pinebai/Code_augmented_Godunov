% 	program velocity

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,3);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\Shu\exact_solution\Shu_v_exact.dat
    xx=Shu_v_exact(:,1);
    yy=Shu_v_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-5.0 5.0 -0.05 3.0]);
   