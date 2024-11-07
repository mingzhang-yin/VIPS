clear; clc;


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
        
        %start updating
        for tries=1:20
            
            %%%%%%%%%%%%%%%
            %generate graph
            %%%%%%%%%%%%%%%
            [A,true_idx,PM,Z]=create_block_model(2000,1,[p_true q_true;q_true p_true],[prop 1-prop]);
            N = size(A,1);
            
            r = randperm(N); % permute row numbers  
            %r = 1:N;
            r_c = sortrows([1:N;r]',2);  r2 = r_c(:,1)';
            A1 = A(r,:); A2 = A1(:,r);

            Z2 = Z(r,:);
            v2 = Z2(:,1) - Z2(:,2);
            
            %%%%%%%%%%%%%%%
            %Initialization
            %%%%%%%%%%%%%%%            
            t = 0.5*log( (p_test*(1-q_test)) / (q_test*(1-p_test)) );
            lambda = 1/(2*t) * log((1-q_test)/(1-p_test)); 
            pi_init = binornd(1,0.5,N,1);
            prop_est = prop;
            
                   
            [pi, P1, P2] = mf(A2,pi_init,S,lambda,t,prop_est,est_flag, v2);         
            ind = (pi>0.5)'+1;
            ind = ind(:,r2);
            ind_onehot = full(ind2vec(ind,2))'; 
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

p_true = 0.2;
q_true = 0.1;
test_p = 0.05:0.03:0.4;
test_q = 0.05:0.03:0.4;

id_p = ((p_true+0.03)-0.05)/0.03;
id_q = ((q_true+0.03)-0.05)/0.03;

[a,b] = size(nmi_r);
for i = 1:a
    for j = i+1:a
        nmi_r(i,j) = -0.1;
    end
end

[a,b] = size(nmi_r);
for i = 1:a
    for j = 1:i-1
        nmi_r(i,j) = nmi_r(i,j)+0.01;
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
title('MFVI','FontSize',20)







