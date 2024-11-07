function r = nmi_perm(Z,ind_onehot)
a = perms([1,2,3]);
n = size(a,1);
nmm = [];
for i = 1:n
    nmm = [nmm, nmi(Z,ind_onehot(:,a(i,:)))];
end

r = max(nmm);

end