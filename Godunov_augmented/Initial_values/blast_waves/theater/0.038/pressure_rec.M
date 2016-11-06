% 	program density

	load d:\Godunov_augmented\show\reconstruction.dat;
    
    x=reconstruction(:,1);
    y=reconstruction(:,6);
    z=reconstruction(:,7);
    
    load d:\Godunov_augmented\show\grid.dat;
    nn=grid(1);
    dd=0.5*(x(2)-x(1));
            
    for i=1: nn
        aa(1)=x(i)-dd;
        aa(2)=x(i);
        bb(1)=y(i)-z(i);
        bb(2)=y(i)-z(i); 
        plot(aa,bb,'r');
        hold on;
        aa(1)=x(i);
        aa(2)=x(i);
        bb(1)=y(i)-z(i);
        bb(2)=y(i)+z(i);
        plot(aa,bb,'r');
        hold on        
        aa(1)=x(i);
        aa(2)=x(i)+dd;
        bb(1)=y(i)+z(i);
        bb(2)=y(i)+z(i);
        plot(aa,bb,'r');
        hold on
    end  
    
    hold on;
    
    plot(x,y,'*b');
    
    hold on;
    
	hold off;
            
    axis([-0.01 1.01 0.0 1.1]);
   