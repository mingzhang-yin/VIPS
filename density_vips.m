clear; clc

prop = 0.3;       

S = 20; 

est_flag = 0
check = false;



p_true = 0.2;
q_true = 0.1;
test_p = 0.05:0.03:0.4;
test_q = 0.05:0.03:0.4;

id_p = ((p_true+0.03)-0.05)/0.03;
id_q = ((q_true+0.03)-0.05)/0.03;


nmi_r = zeros(length(test_p),length(test_q));
nmi_med = zeros(length(test_p),length(test_q));
std_r = zeros(length(test_p),length(test_q));

ii = 1; jj = 1;



for p_test = test_p
    for q_test = test_q
        if p_test <= q_test
            continue
        end
        
        
        NMI = [];

        for tries=1:20
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Generate graph
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    
        [A,true_idx,PM,Z]=create_block_model(2000,1,[p_true q_true;q_true p_true],[prop 1-prop]);
       
        N = size(A,1)/2;
        half = 0.5*ones(N, 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Pairwise partition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        r = randperm(2*N); % permute row numbers
        r_c = sortrows([1:2*N;r]',2);
        r2 = r_c(:,1)';

        A1 = A(r,:);
        A2 = A1(:,r);

        Z2 = Z(r,:);
        zc1 = Z2(1:N,1);
        zc2 = Z2(1:N,2);
        zcp1 = Z2(N+1:2*N,1);
        zcp2 = Z2(N+1:2*N,2);
        v21 = zc1 - zc2; 
        v22 = zcp1 - zcp2;
        v2 = [v21;v22];
        z_star = [zc1; zcp1];

        Azz = A2(1:N,1:N); Ayy = A2(N+1:end,N+1:end);
        Azy = A2(1:N,N+1:end); Ayz = A2(N+1:end,1:N);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Initialize phi, xi, theta10, theta01, theta11, p, q
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        phi_init = binornd(1,0.5,N,1);        
        xi_init = binornd(1,0.5,N,1);
        

        phi_old = phi_init;
        xi_old = xi_init;
        psi10 = zeros(N,1);
        psi01 = zeros(N,1);
        psi11 = zeros(N,1);
        theta10 = normrnd(0,1,N,1)*0;
        theta01 = normrnd(0,1,N,1)*0;
        theta11 = normrnd(0,1,N,1)*0;
        
 
        t = 0.5*log(p_test*(1-q_test)/(q_test*(1-p_test)));
        lambda = 1/(2*t) * log((1-q_test)/(1-p_test));    
        Mzz = 4*t*(Azz - lambda*(ones(N) - eye(N)));
        Myy = 4*t*(Ayy - lambda*(ones(N) - eye(N)));
        Myz = 4*t*(Ayz - lambda*(ones(N) - eye(N)) - diag(diag(Azy)));
        Mzy = 4*t*(Azy - lambda*(ones(N) - eye(N)) - diag(diag(Azy)));
        Szy = 2*t*(diag(Azy) - lambda);
        Syz = 2*t*(diag(Ayz) - lambda);
        
        prop_est = prop;
        


        P1 = [];  %record <u,v1>
        P2 = [];  %record <u,v2>

        %start updating
        for i = 1:S

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update theta10, (phi,xi)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            P1 = [P1, sum([phi_old; xi_old])-N];
            P2 = [P2, [phi_old; xi_old]' * v2];

            theta10 = Mzz*(phi_old - half) + Mzy*(xi_old-half) - Szy + log(prop_est/(1-prop_est));
            [psi10, psi01, psi11, psi00] = update_prob1(theta10,theta01,theta11);
            phi_old = psi10 + psi11;
            xi_old = psi01 + psi11;

            phi = phi_old;
            xi = xi_old; 


            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update theta01, (phi,xi)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            P1 = [P1, sum([phi_old; xi_old])-N];
            P2 = [P2, [phi_old; xi_old]' * v2];

            theta01 = Myy*(xi_old - half) + Myz*(phi_old-half) - Syz + log(prop_est/(1-prop_est));    
            [psi10, psi01, psi11, psi00] = update_prob1(theta10,theta01,theta11);
            phi_old = psi10 + psi11;
            xi_old = psi01 + psi11;

            phi = phi_old;
            xi = xi_old;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update theta11, (phi,xi)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            P1 = [P1, sum([phi_old; xi_old])-N];
            P2 = [P2, [phi_old; xi_old]' * v2];   

            theta11 = (Mzz + Myz)*(phi_old - half) + (Myy + Mzy)*(xi_old - half) + 2*log(prop_est/(1-prop_est));       
            [psi10, psi01, psi11, psi00] = update_prob1(theta10,theta01,theta11);
            phi_old = psi10 + psi11;
            xi_old = psi01 + psi11;  
       
        if i > 3
            [lambda,t,~,~,~] = est_params_pwvb(Azz, Azy, Ayy, Ayz, psi10, psi01, psi11);
            Mzz = 4*t*(Azz - lambda*(ones(N) - eye(N)));
            Myy = 4*t*(Ayy - lambda*(ones(N) - eye(N)));
            Myz = 4*t*(Ayz - lambda*(ones(N) - eye(N)) - diag(diag(Azy)));
            Mzy = 4*t*(Azy - lambda*(ones(N) - eye(N)) - diag(diag(Azy)));
            Szy = 2*t*(diag(Azy) - lambda);
            Syz = 2*t*(diag(Ayz) - lambda);
        end
        
        end
        
    
        phi = phi_old;
        xi = xi_old;    
        label = [phi;xi]';   
        labels(tries,:) = label(:,r2);

        ind = [phi>0.5;xi>0.5]'+1;
        ind = ind(:,r2);
        ind_onehot = full(ind2vec(ind,2))'; %N*2
        NMI = [NMI,nmi(Z,ind_onehot)];
        end
   
    
    nmi_r(ii,jj) =  mean(NMI);
    nmi_med(ii,jj) =  median(NMI);
    std_r(ii,jj) =  std(NMI);
    [p_test,q_test,mean(NMI)]
    
    jj = jj + 1;
    end
    
    ii = ii+1;
    jj = 1;

end




[a,b] = size(nmi_r);
for i = 1:a
    for j = i+1:a
        nmi_r(i,j) = -1;
    end
end

figure;
imagesc(nmi_r');


lowestValue = min(nmi_r(nmi_r(:)>=0))
highestValue = max(nmi_r(:))
cmap = jet(256);
caxis(gca,[lowestValue-2/256, highestValue]);
% Make less than lowest value black:
cmap(1,:)=[211/256,211/256,211/256];
colormap(cmap(1:160,:))
colorbar

set(gca,'YDir','normal') 
set(gca, 'YTick', 1:length(test_q), 'YTickLabel', test_q(1:end)) 
set(gca, 'XTick', 1:length(test_p), 'XTickLabel', test_p(1:end))
hold on
line([id_p,id_p], [0,length(test_q)+1], 'Color', 'k', 'LineStyle','--','LineWidth',2);
line([0,length(test_p)+1],[id_q, id_q], 'Color', 'k', 'LineStyle','--','LineWidth',2);
colorbar()
caxis([0,1])

ax = gca;
ax.FontSize = 14; 
xlabel('$\hat{p}$','Interpreter','latex', 'FontSize',20,'FontWeight','bold')
ylabel('$\hat{q}$','Interpreter','latex',  'FontSize',20,'FontWeight','bold')








