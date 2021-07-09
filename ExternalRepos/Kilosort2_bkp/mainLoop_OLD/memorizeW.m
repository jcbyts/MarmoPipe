function rez = memorizeW(rez, W, dWU,U,mu)

rez.dWU = dWU;
rez.W = gather(W);
rez.U = gather(U);
rez.mu = gather(mu);

rez.Wraw = [];
for n = 1:size(U,2)
    % temporarily use U rather Urot until I have a chance to test it
    rez.Wraw(:,:,n) = rez.mu(n) * sq(rez.U(:,n,:)) * sq(rez.W(:,n,:))';
end
