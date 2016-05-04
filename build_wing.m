function [wing] = build_wing(c,structure,theta)
    %assigns coordinates to the vertices of each skin, spar and bracket
    %component in the wing section

    bl = structure.bl;
    bt = structure.bt;
    sh = structure.sh;
    st = structure.st;
    kt = structure.kt;

    %slanted section length
    X = (-c*cosd(theta)+sqrt(c^2*cosd(theta)^2+3*c^2))/2;
    %wing skins
    wing.skin(1).x = [0 0 c/2 c/2];
    wing.skin(1).y = [sh/2+kt sh/2 sh/2+kt sh/2];
    
    wing.skin(2).x = [c/2 c/2  X*cosd(theta)+c/2 X*cosd(theta)+c/2];
    wing.skin(2).y = [sh/2+kt sh/2 sh/2+kt-X*sind(theta) sh/2-X*sind(theta)];
    
    wing.skin(3).x = [0 0 c/2 c/2];
    wing.skin(3).y = [-sh/2 -sh/2-kt -sh/2 -sh/2-kt];
    
    wing.skin(4).x = [c/2 c/2 X*cosd(theta)+c/2 X*cosd(theta)+c/2];
    wing.skin(4).y = [-sh/2 -sh/2-kt -sh/2-X*sind(theta) -sh/2-kt-X*sind(theta)];
    
    %spars
    wing.spar(1).x = [0 0 st st];
    wing.spar(1).y = [sh/2 -sh/2 sh/2 -sh/2];
    
    wing.spar(2).x = [c/2-st c/2-st c/2 c/2];
    wing.spar(2).y = [sh/2 -sh/2 sh/2 -sh/2];
    
    wing.spar(3).x = [c/2+X*cosd(theta)-st*cosd(theta) c/2+X*cosd(theta)-st*cosd(theta) ...
        c/2+X*cosd(theta) c/2+X*cosd(theta)];
    wing.spar(3).y = [-X*sind(theta)+sh/2+st*sind(theta) -X*sind(theta)-sh/2+st*sind(theta) ...
        -X*sind(theta)+sh/2 -X*sind(theta)-sh/2];
    
    %brackets
    wing.bracket(1).x = [st st st+bt st+bt];
    wing.bracket(1).y = [sh/2 sh/2-bl sh/2 sh/2-bl];
    
    wing.bracket(2).x = [st+bt st+bt st+bt+bl st+bt+bl];
    wing.bracket(2).y = [sh/2 sh/2-bt sh/2 sh/2-bt];
    
    wing.bracket(3).x = [c/2-(st+bt+bl) c/2-(st+bt+bl) c/2-(st+bt) c/2-(st+bt)];
    wing.bracket(3).y = [sh/2 sh/2-bt sh/2 sh/2-bt];
    
    wing.bracket(4).x = [c/2-(st+bt) c/2-(st+bt) c/2-st c/2-st];
    wing.bracket(4).y = [sh/2 sh/2-bl sh/2 sh/2-bl];
    
    wing.bracket(5).x = [c/2 c/2 c/2+bt*cosd(theta) c/2+bt*cosd(theta)];
    wing.bracket(5).y = [sh/2 sh/2-bl sh/2-bt*sind(theta) sh/2-(bt*sind(theta)+bl)];
    
    wing.bracket(6).x = [c/2+bt*cosd(theta) c/2+bt*cosd(theta) c/2+(bt+bl)*cosd(theta) ...
        c/2+(bt+bl)*cosd(theta)];
    wing.bracket(6).y = [sh/2-bt*sind(theta) sh/2-(bt*sind(theta)+bt) ...
        sh/2-(bt+bl)*sind(theta) sh/2-(bt+bl)*sind(theta)-bt];
    
    wing.bracket(7).x = [c/2+(X-bl-bt-st)*cosd(theta) c/2+(X-bl-bt-st)*cosd(theta) ...
        c/2+(X-bt-st)*cosd(theta) c/2+(X-bt-st)*cosd(theta)];
    wing.bracket(7).y = [sh/2-(X-bl-bt-st)*sind(theta) sh/2-(X-bl-bt-st)*sind(theta)-bt ...
        sh/2-(X-bt-st)*sind(theta) sh/2-(X-bt-st)*sind(theta)-bt];
    
    wing.bracket(8).x = [c/2+(X-st-bt)*cosd(theta) c/2+(X-st-bt)*cosd(theta) ...
        c/2+(X-st)*cosd(theta) c/2+(X-st)*cosd(theta)];
    wing.bracket(8).y = [sh/2-(X-st-bt)*sind(theta) sh/2-(X-st-bt)*sind(theta)-bl ...
        sh/2-(X-st)*sind(theta) sh/2-(X-st)*sind(theta)-bl];
    
    wing.bracket(9).x = [st st st+bt st+bt];
    wing.bracket(9).y = [bl-sh/2 -sh/2 bl-sh/2 -sh/2];
    
    wing.bracket(10).x = [st+bt st+bt st+bt+bl st+bt+bl];
    wing.bracket(10).y = [-sh/2+bt -sh/2 -sh/2+bt -sh/2];
    
    wing.bracket(11).x = [c/2-bt-bl-st c/2-st-bt-bl c/2-bt-st c/2-st-bt];
    wing.bracket(11).y = [-sh/2+bt -sh/2 -sh/2+bt -sh/2];
    
    wing.bracket(12).x = [c/2-st-bt c/2-st-bt c/2-st c/2-st];
    wing.bracket(12).y = [-sh/2+bl -sh/2 -sh/2+bl -sh/2];
    
    wing.bracket(13).x = [c/2 c/2 c/2+bt*cosd(theta) c/2+bt*cosd(theta)];
    wing.bracket(13).y = [-sh/2+bl -sh/2 -sh/2+bl-bt*sind(theta) -sh/2-bt*sind(theta)];
    
    wing.bracket(14).x = [c/2+bt*cosd(theta) c/2+bt*cosd(theta) ...
        c/2+bt*cosd(theta)+bl*cosd(theta) c/2+bt*cosd(theta)+bl*cosd(theta)];
    wing.bracket(14).y = [-sh/2+bt-bt*sind(theta) -sh/2-bt*sind(theta) ...
        -sh/2+bt-bt*sind(theta)-bl*sind(theta) -sh/2-bt*sind(theta)-bl*sind(theta)];

    wing.bracket(15).x = [c/2+X*cosd(theta)-(st+bt+bl)*cosd(theta) ...
        c/2+X*cosd(theta)-(st+bt+bl)*cosd(theta) c/2+X*cosd(theta)-(st+bt)*cosd(theta) ...
        c/2+X*cosd(theta)-(st+bt)*cosd(theta)];
    wing.bracket(15).y = [-X*sind(theta)-sh/2+(st+bt+bl)*sind(theta)+bt ...
        -X*sind(theta)-sh/2+(st+bt+bl)*sind(theta) ...
        -X*sind(theta)-sh/2+bt+(st+bt)*sind(theta) ...
        -X*sind(theta)-sh/2+(st+bt)*sind(theta)];
    
    wing.bracket(16).x = [c/2+(X-st-bt)*cosd(theta) c/2+(X-st-bt)*cosd(theta) ...
        c/2+(X-st)*cosd(theta) c/2+(X-st)*cosd(theta)];
    wing.bracket(16).y = [(st+bt-X)*sind(theta)+bl-sh/2 (st+bt-X)*sind(theta)+bl-sh/2-bl ...
        (st-X)*sind(theta)+bl-sh/2 (st-X)*sind(theta)-sh/2];
    
%     hold on;
%     plot(wing.skin(1).x,wing.skin(1).y);
%     plot(wing.skin(2).x,wing.skin(2).y);
%     plot(wing.skin(3).x,wing.skin(3).y);
%     plot(wing.skin(4).x,wing.skin(4).y);
%     
%     plot(wing.spar(1).x,wing.spar(1).y);
%     plot(wing.spar(2).x,wing.spar(2).y);
%     plot(wing.spar(3).x,wing.spar(3).y);
%     
%     plot(wing.bracket(1).x,wing.bracket(1).y);
%     plot(wing.bracket(2).x,wing.bracket(2).y);
%     plot(wing.bracket(3).x,wing.bracket(3).y);
%     plot(wing.bracket(4).x,wing.bracket(4).y);
%     plot(wing.bracket(5).x,wing.bracket(5).y);
%     plot(wing.bracket(6).x,wing.bracket(6).y);
%     plot(wing.bracket(7).x,wing.bracket(7).y);
%     plot(wing.bracket(8).x,wing.bracket(8).y);
%     plot(wing.bracket(9).x,wing.bracket(9).y);
%     plot(wing.bracket(10).x,wing.bracket(10).y);
%     plot(wing.bracket(11).x,wing.bracket(11).y);
%     plot(wing.bracket(12).x,wing.bracket(12).y);
%      plot(wing.bracket(13).x,wing.bracket(13).y);
%     plot(wing.bracket(14).x,wing.bracket(14).y);
%     plot(wing.bracket(15).x,wing.bracket(15).y);
%     plot(wing.bracket(16).x,wing.bracket(16).y);
%     axis equal
%     
end