% 	program density

	load d:\Godunov_augmented\show\solution.dat;
    
    x=solution(:,1);
    y=solution(:,2);
    
	plot(x,y,'Or'); 
    hold on;
    
    load d:\Godunov_augmented\Initial_values\blast_waves\theater\0.038\exact_solution\Blast_d_exact.dat
    xx=Blast_d_exact(:,1);
    yy=Blast_d_exact(:,2);
    plot(xx,yy,'-');
    hold off
            
    axis([-0.0 1.0 0.0 7.0]);
   