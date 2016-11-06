% 	program velocity

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,3);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\show\exact\exact_solution.dat
    xx=exact_solution(:,1);
    yy=exact_solution(:,3);
    plot(xx,yy,'-');
    hold off
            
    axis([-0.1 1.1 0.0 1.2]);
   