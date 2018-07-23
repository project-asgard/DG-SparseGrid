function pde=Maxwell1
%====================================
% Numerical Test 1 for Maxwells' Eq.
%====================================


c0=3e8;
w=pi/2;%2*pi*sqrt(3)/2*c0;
mu=1;
eps=1;

    function E=exact_E(x,y,z,t)
        E(:,1)=-cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z).*cos(w*t);
        E(:,2)= x-x;
        E(:,3)= sin(2*pi*x).*sin(2*pi*y).*cos(2*pi*z).*cos(w*t);
    end

    function E=Ex(x)% without t
        E(:,1)=-cos(2*pi*x);
        E(:,2)= x-x;
        E(:,3)= sin(2*pi*x);
    end
    function E=Ey(y)% without t
        E(:,1)=sin(2*pi*y);
        E(:,2)=y-y;
        E(:,3)=sin(2*pi*y);
    end
    function E=Ez(z)% without t
        E(:,1)=sin(2*pi*z);
        E(:,2)=z-z;
        E(:,3)=cos(2*pi*z);
    end

    function B=exact_B(x,y,z,t)
        B(:,1)=-2*pi/w*sin(2*pi*x).*cos(2*pi*y).*cos(2*pi*z).*sin(w*t);
        B(:,2)= 4*pi/w*cos(2*pi*x).*sin(2*pi*y).*cos(2*pi*z).*sin(w*t);
        B(:,3)=-2*pi/w*cos(2*pi*x).*cos(2*pi*y).*sin(2*pi*z).*sin(w*t);
    end

    function B=Bx(x)
       B(:,1) = -2*pi/w*sin(2*pi*x);
       B(:,2) =  4*pi/w*cos(2*pi*x);
       B(:,3) = -2*pi/w*cos(2*pi*x);
    end
    function B=By(y)
        B(:,1) = cos(2*pi*y);
        B(:,2) = sin(2*pi*y);
        B(:,3) = cos(2*pi*y);
    end
    function B=Bz(z)
        B(:,1) = cos(2*pi*z);
        B(:,2) = cos(2*pi*z);
        B(:,3) = sin(2*pi*z);
    end

    function f=rhs1(x,y,z,t)
        f(:,1)=x-x;
        f(:,2)=y-y;
        f(:,3)=z-z;
    end

    function f=rhs2(x,y,z,t)
        f(:,1)=-eps*cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z).*sin(w*t)*w...
                +12*pi^2*cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z).*sin(w*t)/w/mu;
        f(:,2)= x-x;
        f(:,3)= eps*sin(2*pi*x).*sin(2*pi*y).*cos(2*pi*z).*sin(w*t)*w...
            -12*pi^2*sin(2*pi*x).*sin(2*pi*y).*cos(2*pi*z).*sin(w*t)/w/mu;
    end
% let rhs2=rhs3*sin(w*t)+rhs4*sin(w*t)

    function f=rhs3x(x)
        f(:,1)=-eps*cos(2*pi*x)*w;
        f(:,2)=x-x;
        f(:,3)= eps*sin(2*pi*x)*w;
    end

    function f=rhs3y(y)
        f(:,1)=sin(2*pi*y);
        f(:,2)=y-y;
        f(:,3)= sin(2*pi*y);
    end

    function f=rhs3z(z)
        f(:,1)=sin(2*pi*z);
        f(:,2)=z-z;
        f(:,3)=cos(2*pi*z);
    end

    function f=rhs4x(x)
        f(:,1)= 12*pi^2*cos(2*pi*x)/w/mu;
        f(:,2)=x-x;
        f(:,3)=-12*pi^2*sin(2*pi*x)/w/mu;
    end

    function f=rhs4y(y)
        % sign
        f(:,1)=sin(2*pi*y);
        f(:,2)=y-y;
        f(:,3)=sin(2*pi*y);
    end

    function f=rhs4z(z)
        f(:,1)=sin(2*pi*z);
        f(:,2)=z-z;
        f(:,3)=cos(2*pi*z);
    end

pde = struct('mu',mu,'eps',eps,'w',w,'E',@exact_E,'B',@exact_B,'fE',@rhs1,'fB',@rhs2,...
    'rhs3x',@rhs3x,'rhs3y',@rhs3y,'rhs3z',@rhs3z,'rhs4x',@rhs4x,'rhs4y',@rhs4y,'rhs4z',@rhs4z,...
    'Ex',@Ex,'Ey',@Ey,'Ez',@Ez,'Bx',@Bx,'By',@By,'Bz',@Bz);

end