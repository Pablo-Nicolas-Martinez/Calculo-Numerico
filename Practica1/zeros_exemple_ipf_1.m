% Calcul del zero de x=exp(-x) per iteració simple.
% a I=[exp(-1),log(2)].

clear;
format long e;
format compact;

% Valors inicials.
I=[1/2,log(2)]; a=0.5; tol=1.e-4; itmax=100;

figure(1);
clf(1);

I=[0,1];
% Valors als extrems de I de g.
fprintf('\n');
fprintf('     I(1)=%24.16e      I(2)=%24.16e \n',I(1), I(2));
fprintf('  g(I(1))=%24.16e   g(I(2))=%24.16e \n',g(I(1)), g(I(2)));
% Valor derivada g.
fprintf(' gp(I(1))=%24.16e  gp(I(2))=%24.16e \n',gp(I(1)), gp(I(2)));

I=[1/2,log(2)];
% Valors als extrems de I de g.
fprintf('\n');
fprintf('     I(1)=%24.16e      I(2)=%24.16e \n',I(1), I(2));
fprintf('  g(I(1))=%24.16e   g(I(2))=%24.16e \n',g(I(1)), g(I(2)));
% Valor derivada g.
fprintf(' gp(I(1))=%24.16e  gp(I(2))=%24.16e \n',gp(I(1)), gp(I(2)));


% Calculem primer la solució amb molta precissió.
[xk,res1,it1] = fixed_point_iter(a,1.e-12,itmax,@g);
alpha=xk(end);
fprintf('\n');
fprintf(' Solució calculada amb molt precissió (tol=1.e-12). \n');
fprintf(' Solucio=%24.16e  Nombre iteracions=%3i \n',xk(end), it1);
fprintf(' gp(alpha)=%24.16e  \n',gp(xk(end)));

L=exp(-I(1));
fprintf('\n');
fprintf(" Màxim de gp  l'interval. \n");
fprintf(' L=%24.16e \n',L);


iter=log((1-L)*tol/abs(exp(-a)-a))/log(L);
fprintf('\n');
fprintf(' Predicció nombre iteracions=%3i \n', iter);
fprintf('\n');
iter=ceil(iter);
fprintf(' Predicció nombre iteracions=%3i \n', iter);

% Calculem la solució amb tol.
fprintf('\n');
[xk,res1,it1] = fixed_point_iter(a,tol,itmax,@g);
fprintf(' Solucio=%24.16e  Nombre iteracions=%3i \n',xk(end), it1);

l1=size(xk,2);
semilogy((0:l1-1),abs(xk-alpha),'r-o');
xlabel('k');
ylabel('|x_k-\alpha|');

hold on;
semilogy((0:size(xk,2)-1),tol*ones(1,size(xk,2)),'k-');
hold off;

figure(2);
clf(2);
plot((0:size(xk,2)-1),xk,'b-o');
%plot((0:size(xk,2)-1),xk-alpha,'b-o');
xlabel('k');
ylabel('x_k');
%ylabel('x_k-\alpha');


figure(3);
clf(3);
semilogy((1:size(xk,2)-1),(xk(2:end)-alpha)./(xk(1:end-1)-alpha),'k-o');
xlabel('k');
ylabel('(x_k-\alpha)/(x_{k-1}-\alpha)');

% Valor derivada g a l''alpha
fprintf('\n\n');
fprintf('                         gp(alpha)=%24.16e  \n',gp(alpha));
fprintf(' (xk(end)-alpha)/(xk(end-1)-alpha)=%24.16e  \n', ...
	 (xk(end)-alpha)/(xk(end-1)-alpha));

figure(4);
clf(4);

%semilogy((1:size(xk,2)),res,'r-o');
semilogy((1:size(xk,2)-1),abs(xk(2:end)-xk(1:end-1)),'r-o');
xlabel('k');
ylabel('|x_k-x_{k-1}|');
hold on;
semilogy((1:size(xk,2)),tol*ones(1,size(xk,2)),'k-');
hold off;

function f = g(x)
 f = exp(-x);
end

function f = gp(x)
 f = -exp(-x);
end
