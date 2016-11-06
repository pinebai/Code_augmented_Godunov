% 	program entropy

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,5);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\Shu\exact_solution\Shu_e_exact.dat
    xx=Shu_e_exact(:,1);
    yy=Shu_e_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-5.0 5.0 -0.5 2.6]);
   