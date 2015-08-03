function BM_02_08_missinglinefix(t,r,N)

    rng('shuffle');

    eta = 1.0e-3; % viscosity of water in SI units (Pascal-seconds)
    kB = 1.38e-23; % Boltzmann constant
    T = 293; % Temperature in degrees Kelvin
    tau = .1; % time interval in seconds
    d_t=5.1e-8; % diameter in meters of toehold
    d_r=1.7e-8; % of RNA
    d_c=5.1e-8; % of toehold-RNA complex
    D_t = kB * T / (3 * pi * eta * d_t); %diffusion coefficient
    D_r = kB * T / (3 * pi * eta * d_r);
    D_c = kB * T / (3 * pi * eta * d_c);
    p_t = sqrt(2*D_t*tau); 
    p_r = sqrt(2*D_r*tau);
    p_c = sqrt(2*D_c*tau);
    
    %CONFINEMENT - height can be changed depending on the volume of the solution (rather than the total possible volume of the eppendorf)
    A = 2.5e-2; %binding distance, default at 1e-7, real world value 2.5e-10
    cylinder_radius=4.6e-3; %radius of 1.5mm eppendorf in metres (3.64e-3metres)
    tube_height_max=18e-3; %height of 1.5mm eppendorf (36e-3metres total)
    tube_height_min=-18e-3;
    cone_height=18e-3;

    %GIF stuff
    theta = 0; % changes the viewing angle
    
    %Figure Stuff
    figure()
    %axis([-0.00005 0.00005 -0.00005 0.00005 -0.00005 0.00005]);
    axis([-4.6e-3 4.6e-3 -4.6e-3 4.6e-3  -18e-3 18e-3])
    grid on
    grid MINOR
    set(gcf, 'Position', [100 10 600 600])
    xlabel('X-axis')
    ylabel('Y-axis')
    zlabel('Z-axis')
    hold on

    if t>=r
        c=r;
    else
        c=t;
    end
    
    blank=[0,0,0];
    points={};
    
    for j=1:N
        for k=1:t+r+(2*c)
            points{j,k}=blank;
        end
        for k=t+r+2:2:t+r+(2*c)
            points{j,k}=0;
        end
    end

    for j=1:N
        if j==1
            [coords, startposition] = startpoint();
%             c_joined=zeros(1,c);
%             c_split=zeros(1,c);
            check=zeros(c,3);
        else
            %Toehold
            for i=1:t
                if any(any(check==i))~=1
                    for k=1:3
                        coords(k+3,i)=coords(k,i);
                        coords(k,i)=cumsum(p_t*randn(1))+startposition(k,i);
                    end   

                    if coords(3,i)>=tube_height_max
                        coords(3,i)=tube_height_max;
                    elseif coords(3,i)<=tube_height_min
                        coords(3,i)=tube_height_min;
                    end

                    if coords(3,i)>=0
                        [coords]=checkxy(cylinder_radius, coords, i);
                    else
                        cone_radius=(0.26*coords(3,i))+4.6e-3;
                        [coords]=checkxy(cone_radius, coords, i);
                    end
                end
            end
            %RNA
            for i=t+1:t+r
                if any(any(check==i))~=1
                    for k=1:3
                        coords(k+3,i)=coords(k,i);
                        coords(k,i)=cumsum(p_r*randn(1))+startposition(k,i);
                    end   
                    if coords(3,i)>=tube_height_max
                        coords(3,i)=tube_height_max;
                    end
                    if coords(3,i)<=tube_height_min
                        coords(3,i)=tube_height_min;
                    end
                    if coords(3,i)>=0
                        [coords]=checkxy(cylinder_radius, coords, i);
                    else
                        cone_radius=(0.26*coords(3,i))+4.6e-3;
                        [coords]=checkxy(cone_radius, coords, i);
                    end
                end
            end
            %Complex
            for i=t+r+1:t+r+c
                if check(i-(t+r),1)~=0 && check(i-(t+r),2)~=0
                    for k=1:3
                        coords(k+3,i)=coords(k,i);
                        coords(k,i)=cumsum(p_c*randn(1))+startposition(k,i);
                    end
                    if coords(3,i)>=tube_height_max
                        coords(3,i)=tube_height_max;
                    end
                    if coords(3,i)<=tube_height_min
                        coords(3,i)=tube_height_min;
                    end
                    if coords(3,i)>=0
                        [coords]=checkxy(cylinder_radius, coords, i);
                    else
                        cone_radius=(0.26*coords(3,i))+4.6e-3;
                        [coords]=checkxy(cone_radius, coords, i);
                    end
                end
            end
        end
        for q=1:c
            if check(q,3)==0 && check(q,1)==0 && check(q,2)==0 && j~=1
                jointime=j; %variable to make sure joiner and splitter dont happen in same time step 
                [coords,check, startposition, points] = joiner(coords,check, startposition, points);
            elseif check(q,3)~=0
                check(q,3)=check(q,3)+1;
            end
            if check(q,1)~=0 && check(q,2)~=0 && j~=jointime
                [coords, check, startposition, points] = splitter(q, coords, check, startposition, points);
            end
        end
        for f=1:t+r+c
            if f<t+r+1
                matcoords=[coords(1,f),coords(2,f),coords(3,f)];
                points{j,f}=matcoords;
            end
            if f>=(t+r+1)
                if f==(t+r+1)
                    g=f;
                    h=f+1;
                end
                matcoords=[coords(1,f),coords(2,f),coords(3,f)];
                points{j,g}=matcoords;
                g=g+2;
                if j>1
                    if points{j-1,h}(1,1)~=0
                        if points{j,h}==0 
                            points{j,h}=1;
                        end
                    end
                end
                h=g+1;    
            end
        end
        
        if j==N
            plotter(points);
            finish=1;
        end
    end
    
    %creates start point definitely inside dimensions of eppendorf
    function [coords, startposition] = startpoint()
        coords=zeros(6,(t+r+c)); %each column is a toehold, with six rows, for current xyz and previous xyz
%         startposition=zeros(6,(t+r+c));
        for m=1:t+r
            coords(3,m)=(tube_height_min)+((tube_height_max-tube_height_min)*rand(1));
            if coords(3,m)>=0
                coords(1,m)=(-cylinder_radius)+((cylinder_radius-(-cylinder_radius))*rand(1));
                coords(2,m)=(-cylinder_radius)+((cylinder_radius-(-cylinder_radius))*rand(1));
            end
            if coords(3,m)<0
                cone_radius=(0.26*coords(3,m))+4.6e-3;
                coords(1,m)=(-cone_radius)+((cone_radius-(-cone_radius))*rand(1));
                coords(2,m)=(-cone_radius)+((cone_radius-(-cone_radius))*rand(1));
            end
        end
        startposition=coords;
    end
       
    %checks the coordinates are within the boundaries of the eppendorf 
    
    function [coords]=checkxy(radius,coords,i)
        if (coords(1,i)^2)+(coords(2,i)^2)>=(radius^2)
            grad=abs(coords(2,i)/coords(1,i));
            if coords(1,i)<0
               coords(1,i)=-(((radius^2)/((grad^2)+1))^0.5); 
            else
               coords(1,i)=(((radius^2)/((grad^2)+1))^0.5);
            end
            if coords(2,i)<0
               coords(2,i)=-(grad*(((radius^2)/((grad^2)+1))^0.5));
            else
               coords(2,i)=grad*(((radius^2)/((grad^2)+1))^0.5);
            end
        end
    end

    function [coords, check, startposition, points] = joiner(coords, check, startposition, points)
        colshift=0;
        joinprob = randi([1 10],1);
        if joinprob >=3
            for n=1:t
                for m=t+1:r+t
                    if any(any(check==n))==1 || any(any(check==m))==1
                        continue
                    else
       %                 if ((((tx(j,k)-rx(j,m))^2)+((ty(j,k)-ry(j,m))^2)+((tz(j,k)-rz(j,m))^2))<=(A^2) || (check_r(1,m)~=0 && check_t(1,k)~=0)) && (j~=1) && delay(1,n)==0    
                        %dist2 on sub matrices to find distances between for joining 
                        if (((coords(1,n)-coords(1,m))^2)+((coords(2,n)-coords(2,m))^2)+((coords(3,n)-coords(3,m))^2))<=(A^2)
                            for p=1:c    
                                if check(p,1)~=0 && check(p,2)~=0
                                    continue
                                else
                                    if check(p,1)==0 && check(p,2)==0
                                        coords(1,t+r+p)=(coords(1,n)+coords(1,m))/2;
                                        coords(2,t+r+p)=(coords(2,n)+coords(2,m))/2;
                                        coords(3,t+r+p)=(coords(3,n)+coords(3,m))/2;
                                        startposition(1,t+r+p)=coords(1,t+r+p);
                                        startposition(2,t+r+p)=coords(2,t+r+p);
                                        startposition(3,t+r+p)=coords(3,t+r+p);
                                        check(p,1)=n;
                                        check(p,2)=m;
                                        if colshift==0;
                                            colshift=1;
                                        end
                                        points{j,t+r+p+colshift}=[1,n,m];
                                        colshift=colshift+1;

                                    end
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    function [coords, check, startposition, points] = splitter(q, coords, check, startposition, points)
        split = randi([1 10],1);
        if split<=4
            if check(q,1)~=0 && check(q,2)~=0
                toehold=check(q,1);
                rna=check(q,2);
                startposition(1,toehold)=coords(1,t+r+q);
                startposition(2,toehold)=coords(2,t+r+q);
                startposition(3,toehold)=coords(3,t+r+q);
                startposition(1,rna)=coords(1,t+r+q);
                startposition(2,rna)=coords(2,t+r+q);
                startposition(3,rna)=coords(3,t+r+q);
                
                points{j,t+r+(2*q)}=[0,toehold,rna];
                check(q,1)=0;
                check(q,2)=0;
                check(q,3)=-5;
            end
        end
    end

    function [] = plotter(points)
        for row=1:N
            if row==1
                %startpoints
                for col=1:t+r
                   plot3(points{row,col}(1,1), points{row,col}(1,2), points{row,col}(1,3),'kx') 
                end
            else
                for altcol=t+r+2:2:t+r+(2*c)
                    if points{row,altcol}(1,1)==1 && size(points{row,altcol},2)==3 && row>1
                        if points{row-1,altcol}(1,1)==0 %definitely joinpoint
                            %plot blue using row-1 to row to joinpoint column, column selected from points{row,col}(1,2)) 
                            plot3([points{row-1,(points{row,altcol}(1,2))}(1,1), points{row,altcol-1}(1,1)], [points{row-1,(points{row,altcol}(1,2))}(1,2), points{row,altcol-1}(1,2)], [points{row-1,(points{row,altcol}(1,2))}(1,3), points{row,altcol-1}(1,3)], 'b')
                            %plot red using row-1 to row to joinpoint column, column selected from points{row,col}(1,2)
                            plot3([points{row-1,(points{row,altcol}(1,3))}(1,1), points{row,altcol-1}(1,1)], [points{row-1,(points{row,altcol}(1,3))}(1,2), points{row,altcol-1}(1,2)], [points{row-1,(points{row,altcol}(1,3))}(1,3), points{row,altcol-1}(1,3)], 'r')
                            %plot circle at joinpoint
                            plot3(points{row,altcol-1}(1,1), points{row,altcol-1}(1,2), points{row,altcol-1}(1,3), 'ko')
                        end
                        if points{row-1,altcol}(1,1)==1 %defintely splitpoint
                            %do nothing
                        end
                    elseif (points{row,altcol}(1,1)==1 && size(points{row,altcol},2)==1) || (points{row,altcol}(1,1)==0 && size(points{row,altcol},2)==3)
                        %plot green from row-1 to row
                        plot3([points{row-1,altcol-1}(1,1), points{row, altcol-1}(1,1)], [points{row-1,altcol-1}(1,2), points{row, altcol-1}(1,2)], [points{row-1,altcol-1}(1,3), points{row, altcol-1}(1,3)],'g')
                    end
                    if points{row,altcol}(1,1)==0 %currently split
                        if size(points{row-1,altcol},2)==3 && points{row-1,altcol}(1,1)==0 %prevents incorrect plotting in the case of join/split in consecutive timesteps
                            %plot from split column to blue column
                            plot3([points{row-1,altcol-1}(1,1), points{row,(points{row-1,altcol}(1,2))}(1,1)],[points{row-1,altcol-1}(1,2), points{row,(points{row-1,altcol}(1,2))}(1,2)],[points{row-1,altcol-1}(1,3), points{row,(points{row-1,altcol}(1,2))}(1,3)],'b')
                            %plot from split column to red column
                            plot3([points{row-1,altcol-1}(1,1), points{row,(points{row-1,altcol}(1,3))}(1,1)],[points{row-1,altcol-1}(1,2), points{row,(points{row-1,altcol}(1,3))}(1,2)],[points{row-1,altcol-1}(1,3), points{row,(points{row-1,altcol}(1,3))}(1,3)],'r')
                            %plot star at splitpoint
                            plot3(points{row-1,altcol-1}(1,1), points{row-1,altcol-1}(1,2), points{row-1,altcol-1}(1,3), 'k*')
                            %plot 
                            continue %might not be necessary
                        end
                    end
                    for col=1:t+r
                        if points{row,col}~=points{row-1,col}
                           if size(points{row,altcol},2)~=3
                               %determine column in t or r for colour of line
                               if col<=t %plot blue for toeholds
                                   if ((size(points{row-1,altcol},2)==3 && points{row-1,altcol}(1,1)==0) && points{row-1,altcol}(1,2)==col)
                                   %if just split, don't plot from row-1 to row in t column since plotting from split point has just occurred
                                       continue
                                   elseif points{row,altcol}(1,1)==0 && size(points{row-1,altcol},2)~=3
                                       plot3([points{row-1,col}(1,1), points{row,col}(1,1)],[points{row-1,col}(1,2), points{row,col}(1,2)],[points{row-1,col}(1,3), points{row,col}(1,3)],'b')
                                   end
                               end
                               if col>t %plot red for rna(trigger)
                                   if ((size(points{row-1,altcol},2)==3 && points{row-1,altcol}(1,1)==0) && points{row-1,altcol}(1,3)==col)
                                   %if just split, don't plot from row-1 to row in r column since plotting from split point has just occurred
                                       continue
                                   elseif points{row,altcol}(1,1)==0 && size(points{row-1,altcol},2)~=3
                                       plot3([points{row-1,col}(1,1), points{row,col}(1,1)],[points{row-1,col}(1,2), points{row,col}(1,2)],[points{row-1,col}(1,3), points{row,col}(1,3)],'r')
                                   end
                               end
                           end
                        end
                    end
                end
            end
        end
    end  
    function jiff(row) 
       change = 360/N; % the size of the angle change 
       % gif utilities
       set(gcf,'color','w'); % set figure background to white
       drawnow;
       frame = getframe(gcf);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       outfile = '17_07_v1.gif';
       
       % adjusting the viewing the angle
       view(theta,45);
       theta = theta + change;

       % On the first loop, create the file. In subsequent loops, append.
       if row==2
          imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
       else
          imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
       end
    end 
end