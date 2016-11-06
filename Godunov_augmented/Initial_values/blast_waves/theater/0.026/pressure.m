% 	program pressure

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,4);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\blast_waves\theater\0.026\exact_solution\Blast_p_exact.dat
    xx=Blast_p_exact(:,1);
    yy=Blast_p_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-0.0 1.0 -2.0 280.1]);
   