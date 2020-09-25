%glueing them together!
%set params and data:
%Single point continuations
nu = 1;
s = 0;
incplot = 1;
d = 1; 
inttolspec = 1; int_tol = 1e-4;
n = 100; 
xl0 = .3; 
xu0 = .5; 
h0 = s*linspace(0,1,100)+ones(1,100);

regime1_v3_noparams;
post_processing_no_plots;
h_reg1 = h;
t_reg1 = t;
x_reg1 = x;
xu_reg1 = xu; xl_reg1 = xl;

xcu_reg1 = xcu; xcl_reg1 = xcl;
hl_reg1 = hl; hu_reg1 = hu;


if sol.ie == 2 %goes to regime 2;
    prevrun =1;
    int_tol = 1e-7;
    regime2_v3;
    post_processing_no_plots;
    x_reg2 = x;
    xcu_reg2 = xcu; xcl_reg2 = xcl;
    xl_reg2 = xl; xu_reg2 = xu;
    hl_reg2 = hl; hu_reg2 = hu;
    h_reg2 = h;
    t_reg2 = t;

    if sol.ie == 2  %goes to regime 3;
        int_tol = 1e-8;
        regime3_v3;
        post_processing_no_plots;
        h_reg3 = h;
        t_reg3 = t;
        xcu_reg3 = xcu; xcl_reg3 = xcl;
        hl_reg3 = hl; hu_reg3 = hu;
        x_reg3 = x;
        xu_reg3 = xu; xl_reg3 = xl;
        
        
        
        h = [h_reg1; h_reg2; h_reg3]; t = [t_reg1 t_reg2 t_reg3];
        x = [x_reg1; x_reg2; x_reg3];
        xcu = [xcu_reg1; xcu_reg2; xcu_reg3]; xcl = [xcl_reg1; xcl_reg2; xcl_reg3];
        hl = [hl_reg1; hl_reg2; hl_reg3]; hu = [hu_reg1; hu_reg2; hu_reg3];
        xl = [xl_reg1; xl_reg2; xl_reg3]; xu = [xu_reg1; xu_reg2; xu_reg3];
    else
        h = [h_reg1; h_reg2]; t = [t_reg1 t_reg2];
        x = [x_reg1; x_reg2];
        xcu = [xcu_reg1; xcu_reg2]; xcl = [xcl_reg1; xcl_reg2];
        hl = [hl_reg1; hl_reg2]; hu = [hu_reg1; hu_reg2];
        xl = [xl_reg1; xl_reg2]; xu = [xu_reg1; xu_reg2];
    end
    
else
h = h_reg1; t = [t_reg1];
xcu = [xcu_reg1]; xcl = [xcl_reg1];
hl = [hl_reg1]; hu = [hu_reg1];
x = x_reg1;
xl = xl_reg1;xu = xu_reg1;
end


%now plot:
if incplot == 1
for j1 = 1:length(t)-1
            
        figure(3); clf; hold on;
        plot(h(j1,:), x(j1,:),'b',hu(j1,:), xcu(j1,:), 'b', hl(j1,:), xcl(j1,:), 'b',...
            linspace(0, h(j1,1), 100),x(j1,1)*ones(1,100), 'k', linspace(0, h(j1, end), 100), x(j1, end)*ones(1,100), 'k'); 
        %and the negative
        plot(-h(j1,:), x(j1,:),'b',-hu(j1,:), xcu(j1,:), 'b', -hl(j1,:), xcl(j1,:), 'b',...
        -linspace(0, h(j1,1), 100),x(j1,1)*ones(1,100), 'k', -linspace(0, h(j1, end), 100), x(j1, end)*ones(1,100), 'k'); 
        %and the top?
        plot([0,0],[xcu(j1, end) 1],'b')
    axis([-2,2,0,1]); title(['t =' num2str(t(j1))]); ylabel('x');
    pausetime = 0.0004;
    pause(pausetime) %note, time steps not uniform, don't make videos from this
end
end





