var     pih             ${\pi_h}$               (long_name='Domestic inflation')
        %x               $x$                     (long_name='Output gap')
        y               $y$                     (long_name='Output')
        %ynat            ${\bar y}$              (long_name='Natural output')
        %rnat            ${\bar r}$              (long_name='Natural interest rate')
        r               $r$                     (long_name='Nominal interest rate')
        s               $s$                     (long_name='Terms of trade')
        pi              ${\pi}$                 (long_name='CPI inflation')
        %piw             ${\piw}$                (long_name='Wage inflation')
        %p               $p$                     (long_name='CPI level')
        %ph          	${p_h}$                 (long_name='Domestic price level')
        e               $e$                     (long_name='Exchange rate')
        ystar           ${y^*}$                 (long_name='World output')
        pistar          ${\pi^{*}}$             (long_name='World inflation')
        n               ${n}$                   (long_name='Employment')
        %nx              ${nx}$                  (long_name='Net Exports')
        %real_wage 		${w-p}$                 (long_name='Real Wage')
        a               $a$                     (long_name='Risk aversion')
        c               $c$                     (long_name='Domestic consumption')
        %deprec_rate 	$\Delta e_t$            (long_name='Nominal depr. rate')
        m               $m$                     (long_name='Monetary policy')
        %w
        u
        %w_nat
        %w_gap
        %l
        %mu_w
        xhat
        %h
        %ghat
        %ch
        %cf
        %pf
    ;

varexo  eps_star        ${\varepsilon^{*}}$     (long_name='World output shock')
        eps_a           ${\varepsilon^{a}}$     (long_name='World output shock')
        eps_mp           ${\varepsilon^{a}}$    (long_name='Monetary policy shock')
    ;

parameters 
        sigma           $\sigma$                (long_name='risk aversion')
        eta             $\eta$                  (long_name='Substitution home foreign')
        gamma           $\gamma$                (long_name='Substitution between foreign')
        phi             $\varphi$               (long_name='Inverse Frisch elasticity')
        epsilon 		$\varepsilon$           (long_name='Elasticit of substitution')
        theta           $\theta$                (long_name='Calvo parameter')
        beta            $\beta$                 (long_name='discount factor')
        alpha           $\alpha$                (long_name='openness')
        phi_pi          $\phi_\pi$              (long_name='Feedback Taylor rule inflation')
        phi_x           $\phi_x$                (long_name='Feedback Taylor rule output gap')
        rhoa            $\rho_a$                (long_name='autocorrelation TFP')
        rhoy            $\rho_y$                (long_name='autocorrelation foreign output')
        rhom            $\rho_m$                (long_name='autocorrelation monetary policy')
        gam             $\gam$                  (long_name='Hysteresis parameter')
        theta_w         $\theta_w$              (long_name='wage rigidity_Calvo paramater')
        epsilon_w 		$\varepsilon_w$         (long_name='Elasticity of substitution wage')
        lambda_n
        psi_ya
        psi_wa
        lambda_p
        lambda_w
        u_n
        x_n
        delt
        g
        B
        alp
        M
        net0
        net1
        h0
        hL
        hF
        k0
        kL
        kF
        ga
        lam
        psi_p
        cap
        t
        S
        gs
        Psis
        Psia
        h0s
        hLs
        hFs
        omega
        k0s
        kLs
        kFs
        Omega
        
        
    ;

% set deep parameters
sigma = 1;
eta = 1 ;
gamma = 1;
phi = 3;
epsilon = 6;
epsilon_w=4.5;
theta = 0.75;
theta_w=1/24;
beta  = 0.99;
alpha = 0.8;
phi_pi = 1.5;
phi_x = 0.5;
rhoa = 0.9; %use value used for Figure 1, reset later                                                                
rhoy = 0.9;  
rhom = 0.9;
gam=0.9;
alp = 1;

u_n=0.3;
x_n = 0.67;

B=0.12;
ga=0.5;

S=0.5;

lambda_n=(1-theta_w)/(theta_w*epsilon_w);
lambda_w=(1-theta_w)*(1-beta*theta_w)/(theta_w*(1+epsilon_w*phi));

psi_ya=(1+phi)/(sigma*(1-alpha)+phi+alpha);
psi_wa=(1-alpha*psi_ya)/(1-alpha);
lambda_p=((1-theta)*(1-beta*theta)/theta)*((1-alpha)/(1-alpha+alpha*epsilon));

M=epsilon/(epsilon-1);
lam=((1-beta*theta)*(1-theta))/theta;
delt = u_n*x_n/((1-u_n)*(1-x_n));
psi_p = 1-cap*g*M;
t=1-beta*(1-delt);
cap = t;

omega = sigma*gamma+(1-alpha)*(sigma*eta-1);
Omega=omega-1;
g = B*x_n^alp;
gs=g*S^(alpha);
Psia=1-(1-beta*(1-delt))*M*g;
Psis=alpha*(Psia+(1-beta*(1-delt))*M*gs);
net0 = (1-g*(1+alp))/(1-delt*g); 
net1 = g*(1-delt)*(1+alp*(1-x_n))/(1-delt*g);
h0 = (alp*g*M/delt)*(1+beta*(1-delt)*(1-delt)*(1-x_n))+beta*(1-delt)*g*M*(net1-net0);
hL = -(alp*g*M/delt)*(1-delt)*(1-x_n)-beta*(1-delt)*g*M*net1; 
hF = -beta*(1-delt)*g*M*((alp/delt)-net0);

h0s=(Psis*sigma)*(1-net0)/(alpha*omega)+(alp*M*gs)/delt+gs*M*beta*(1-delt)*(net1-net0+alp*(1-delt)*(1-x_n)/delt);
hLs=-Psis*sigma*net1/(alpha*omega)-alp*M*gs*(1-delt)*(1-x_n)/delt-gs*M*beta*(1-delt)*net1;
hFs=-gs*M*beta*(1-delt)*((alp/delt)-net0);

k0 = lam*h0/(1-u_n); 
kL = -lam*hL/(1-u_n);
kF = -lam*hF/(1-u_n);

k0s = lam*h0s/(1-u_n); 
kLs = -lam*hLs/(1-u_n);
kFs = -lam*hFs/(1-u_n);

%kLs=-0.6;
%kFs=0.3;

model;%(linear);

//define parameter dependencies
//steady state real interest rate, defined below equation (11)
#rho  = beta^(-1)-1;
//defined below equation (27)
%#omega = sigma*gamma+(1-alpha)*(sigma*eta-1);
//defined below equation (29)
#sigma_a =sigma/((1-alpha)+alpha*omega);

#Theta=(sigma*gamma-1)+(1-alpha)*(sigma*eta-1);
//defined below equation (32)
#lambda = (1-(beta*theta))*(1-theta)/theta;
//defined below equation (36)
#kappa_a =lambda*(sigma_a+phi);
//defined below equation (35)
#Gamma = (1+phi)/(sigma_a+phi);
#Psi = -Theta*sigma_a/(sigma_a+phi);




%From Gali 2010_Hysteresis

%[name='Wage Phillips curve']
%piw=beta*piw(+1)+(1-gam)*lambda_n*n+gam*lambda_n*(n-n(-1));
%piw=beta*piw(+1)-lambda_w*phi*(u);
%piw=w-w(-1);

%From Blanchard and Gali

%[name='Hiring']
%h=n-n(-1)+delt*n(-1);

%[name='Hiring cost']
%ghat=a+B*alp*xhat;

[name='Labour market tightness']
delt*xhat = n - (1-delt)*(1-x_n)*n(-1); 

%[name='Constraint']
%c=a+((1-g)/(1-delt*g))*n+(g*(1-delt)/(1-delt*g))*n(-1)- (alp*g/(1-delt*g))*delt*xhat;

%[name='Natural level of wage']
%w_nat=psi_wa*a;

%[name='Wage gap']
%w_gap=w_gap(-1)+piw-pih-(w_nat-w_nat(-1));

[name='Unemployment']
%u=l-n;
u=-(1-u_n)*n;
%phi*u=w_gap-(sigma+phi/(1-alpha))*x;

%[name='Wage markup']
%mu_w=phi*u;

%[name='Total comsumption']
%c=((1-alpha)/eta)*ch+(alpha/eta)*cf;

%[name='Home consumption']
%ch=(1-alpha)*(-eta*(ph-p)+c);

%[name='Foreign consumption']
%cf=alpha*(-eta*(pf-p)+c);

%[name='CPI']
%p=(1-alpha)*ph+alpha*(pf);

%[name='Equation (37), IS Curve']
%x    = x(+1) - sigma_a^(-1)*(r - pih(+1) - rnat) ;  
y    = y(+1) - sigma_a^(-1)*(r - pih(+1))+alpha*Omega*(ystar(+1)-ystar);% - rnat) ;                             

[name='Equation (36), Philips Curve']
%pih  = beta * pih(+1)+ kappa_a*x+lambda_p*w_gap;                                                
%pih = beta*pih(+1)-k0*u+kL*u(-1)+kF*u(+1)-lam*Psia*ga*a;
pih = beta*pih(+1)-k0s*u+kLs*u(-1)+kFs*u(+1)-lam*Psia*ga*a;

%[name='Equation below (37)']
%rnat = -sigma_a*Gamma*(1-rhoa)*a + alpha*sigma_a*(Theta+Psi)*(ystar(+1)-ystar);

%[name='Equation (35), definition natural level of output']
%ynat = Gamma*a + alpha*Psi*ystar;                        
                         
%[name='Equation above (35), definition output gap']
%x    = y - ynat;                                   
                            
[name='Equation (29)']
y = ystar + sigma_a^(-1)*s;

[name='Equation (14)']
pi   = pih + alpha*(s-s(-1));

[name='Equation 15 (first difference)']
s    = s(-1) + e - e(-1) + pistar - pih;

[name='Constant world inflation, see p.724 (Given constant world prices)'] 
pistar = 0;

[name='Equation (22), employment']
y = a + n;

%[name='Equation (31), net exports']
%nx = alpha*(omega/sigma-1)*s;

[name='Equation (27), defines consumption']
y = c+alpha*omega/sigma*s;

%[name='Above equation (11), defines real wage']
%w-pi = sigma*c+phi*n;
%real_wage = w-pih;

[name='Equations on p. 723, stochastic processes']
a    = rhoa*a(-1) + eps_a;
ystar= rhoy*ystar(-1) + eps_star;
m    = rhom*m(-1) + eps_mp;

[name='Equations on page'] 
r = phi_pi*pih+phi_x*y+m; // domestic inflation-based Taylor rule (DITR)
  
%[name='definition consumer price level']
%pi   = p - p(-1);

%[name='definition domestic price level']
%pih  = ph - ph(-1);

%[name='definition nominal depreciation rate of exchange rate']
%deprec_rate=e-e(-1);

end;


steady;%(solve_algo=2);
check;

shocks;
var eps_a = 1; //unit shock
%var eps_star = 1;
var eps_mp = 1;
end;

//generate LaTeX-files with equations and parameterization
write_latex_dynamic_model;
write_latex_parameter_table;
write_latex_definitions;

stoch_simul(TeX,order=1,irf=40,irf_plot_threshold=0);

collect_latex_files;

//Uncomment the following lines to generate a PDF file using PDFLaTeX (if installed)
% if system(['pdflatex -halt-on-error ' M_.fname '_TeX_binder.TeX'])
%     error('TeX-File did not compile.')
% end
