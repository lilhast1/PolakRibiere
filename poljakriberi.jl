using Calculus;
using LinearAlgebra
import LinearAlgebra: norm
import Calculus: gradient
import Calculus: hessian
#Numerical Recipies ideje
# 1. Ogradnjivanje minimuma koristeci proceduru Bracketmethod na 491
# 2. Ogradjivanje sa Golden Pretragom
# 3. Linearna pretraga 508
# 4. Sve o poljaku 519

#dodati tex file
#do nedelje bi bilo dobro to rijesiti

function sgn(x) # julia nema built in signum
    return x > 0. ? 1. : -1.;
end

function bracket(a, b, f)
    glimit = 100.;
    phi = (1 + sqrt(5)) / 2.;
    
    ax = a; bx = b;
    fa = f(ax); fb = f(bx);
    if (fb > fa)
        ax, bx = bx, ax;
        fa, fb = fb, fa;
    end
    cx = bx + phi * (bx - ax);
    fc = f(cx);
    while (fb > fc)
        r = (bx - cx) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2. * max(abs(q-r),1e-18) * sgn(q-r));
        ulim = bx + glimit * (cx - bx);
        if ((bx - u) * (u - cx) > 0.)
            fu = f(u);
            if (fu < fc)
                return bx, u, cx;
            elseif (fu > fb)
                return ax, bx, u;
            end
            u = cx + phi * (cx - bx);
            fu = f(u);
        elseif ((cx - u) * (u - ulim) > 0.)
            fu = f(u);
            if (fu < fc)
                bx = cx; fb = fc;
                cx = u; fc = fu;
                u = u + phi * (u - cx);
                fu = f(u);
            end
        elseif ((ulim - cx) * (u - ulim) > 0.)
            u = ulim; fu = f(u);
        else 
            u = cx + phi * (cx - bx);
            fu = f(u);
        end
        ax = bx; fa = fb;
        bx = cx; fb = fc;
        cx = u; fc = fu;
    end
    return ax,bx,cx;
end

function minimize(f, df, ax, bx, cx)
    imax = 200;
    eps = 1e-8;
    
    d = 0.; e = 0.;
    a = min(ax, cx);
    b = max(ax, cx);
    x = bx; v = bx; w = bx;
    fx = f(x); fv = fx; fw = fx;
    dx = df(x); dv = dx; dw = dx;
    
    for i = 1:imax
        xm = (a + b) / 2;
        tol1 = eps * abs(x); tol2 = 2 * tol1;
        if (abs(x - xm) <= tol2 - (b - a) / 2)
            return x, fx;
        end
        if (abs(e) > tol1)
            d1 = 2 * (b - a); d2 = d1;
            if (abs(dw - dx) > eps)
                d1 = (w - x) * dx / (dx - dw); #sekanta
            end
            if (abs(dv - dx) > eps)
                d2 = (v - x) * dx / (dx - dv); #sekanta sa druge strane
            end
            u1 = x + d1;
            u2 = x + d2;
            ok1 = (a - u1) * (u1-b) > 0. && dx * d1 <= 0.;
            ok2 = (a - u2) * (u2-b) > 0. && dx * d2 <= 0.;
            olde = e; e = d;
            if (ok1 || ok2)
                if (ok1 && ok2)
                    d = abs(d1) < abs(d2) ? d1 : d2;
                elseif ok1 
                    d = d1;
                else 
                    d = d2;
                end
                if abs(d) <= abs(olde / 2);
                    u = x + d;
                    if (u - a < tol2 || (b - u) < tol2)
                        d = abs(tol1) * sgn(xm-x);
                    end
                end
            else 
                e = dx > 0. ? a - x : b - x;
                d = e / 2;
            end
        else
            e = dx > 0. ? a - x : b - x;
            d = e / 2;
        end
        if abs(d) > tol1
            u = d + x;
            fu = f(u);
        else 
            u = x + abs(tol1) * sgn(d);
            fu = f(u);
            if (fu > fx)
                return x, fx;
            end
        end
        du = df(u);
        if fu < fx 
            if (u > x)
                a = x;
            else 
                b = x;
            end
            v = w; fv = fw; dv = dw;
            w = x; fw = fx; dw = dx;
            x = u; fx = fu; dx = du;
        else 
            if (u < x)
                a = u;
            else 
                b = u;
            end
            if (fu < fw || abs(x - w) < eps)
                v = w; fv = fw; dv = dw;
                w = u; fw = fu; dw = du;
            elseif (fu < fv || abs(x - v) < eps || abs(v - w) < eps )
                v = u; fv = fu; dv = du;
            end
        end
    end
    return x, fx;
end
x, fx = minimize(x->-exp(-x*x), x->2*x*exp(-x*x), -1, 0.5, 1);
x, fx = minimize(x->x*x, x->2*x, -1,0.6,1);
x, fx = minimize(x->x^3-3*x,x->3*x*x-3,0,0.5,2);
x

a,b,c = bracket(-1, -0.5, x->-exp(-x*x));


function dbrent(a, b, f, df)
    ax = 0.; bx = 0.; cx = 0.;
    # prvo ogradjivanje
    (ax, bx, cx) = bracket(a, b, f);
    # minimizacija unutar ograde
    return minimize(f, df, a, b, c);
end    


function linefactory(f, p, d)
    function line(x)
        return f(p + d * x);
    end
    return line;
end

function direcderivfactoy(f,p,d)
    function direcderiv(x)
        return dot(d, gradient(f, p + x * d));
    end
    return direcderiv;
end

function linmin(f, p, d)
    line = linefactory(f, p, d); # formiramo krivu u pravcu d koju cemo minimizirati
    direcderiv = direcderivfactoy(f,p,d); # formiramo direkcioni izvod
    #pocetno nagadjanje
    a = -1.;
    b = 1.; 
    xmin, fmin = dbrent(a, b, line, direcderiv);  # ova linija bira minimizacijsku shemu po pravcu. 
    # pravljenje koraka
    d = d .* xmin;
    p = p .+ d; 
    return fmin, p, d;
end

function PolakRibiere(f, p)
    imax = 500;
    eps = 1e-10;
    #init proc
    fp = f(p);
    r = -gradient(f, p);
    g = r;
    h = r;
    for i = 1:imax
        fl, p, r = linmin(f, p, r);
        if abs(fl - fp) <= eps * eps 
            return p;
        end
        fp = fl;
        r = -gradient(f, p);
        test = maximum(x->x , abs.(r) .* max.(p, 1) / max(fp, 1));
        
        # if test < eps * eps
        #     return p;
        # end
        gg = dot(g, g);
        dgg = dot(r, r) - dot(g, r);
        if gg == 0.0
            return p;
        end
        b = dgg / g;
        g = r;
        h = g .+ b * h;
        r = h;
    end
    return p;
end



#testiranja

rosenbrock = x->( (1-x[1])^2 + 100* (x[2] - x[1]^2)^2 );
function rastrigin(x) # 5dim, A = 10 min u x=0
    s = 50.;
    for i = 1:5
        s = s + x[i]^2 - 10 * cos(2 * pi * x[i]);
    end
    return s;
end
function ackley(x) #min u 0,0
    return -20 * exp(-0.2 * sqrt((x[1]^2 + x[2]^2) / 2.)) - 
            exp(0.5 * (cos(2 * pi * x[1]) + cos(2 * pi * x[2]))) + exp(1) + 20.;
end
function beale(x) # gloablni min u [3,0.5]
    return (1.5 - x[1] + x[1] * x[2])^2 + (2.25 - x[1] + x[1] * x[2]^2)^2
            + (2.625 - x[1] + x[1] * x[2]^3)^2;
end
#jos dobrih funkcija na https://en.wikipedia.org/wiki/Test_functions_for_optimization
PolakRibiere(x->-exp(-(x[1]*x[1]+x[2]*x[2])), [2.,1.])
PolakRibiere(rosenbrock, [0.5, 0.])
PolakRibiere(x->x[1]^2+x[2]^2, [3.,1.])
PolakRibiere(rastrigin, [1., 1., 1., 1., 1.])
PolakRibiere(ackley, [1.,1.]) # veoma veoma osjetljivo recimo 2,1 eksplodira
PolakRibiere(beale, [2.5,0.2]) # ova je brutalna
PolakRibiere(x->x[1]^4 + x[2]^4 - 4 * x[1] * x[2], [0.5,-0.5]) # 1, 1 ili -1, -1 dok je 0,0 sedlo
PolakRibiere(x->(x[1]^2 + x[2]) * sqrt(exp(x[2])), [0., 0.]) #-2/e je min u 0, -2
PolakRibiere(x->x[1] + x[2]^2/(4*x[1]) + x[3]^2/x[2] + 2/ x[3], [.4, 1.1, 1.1]) # u 0.5,1,1

