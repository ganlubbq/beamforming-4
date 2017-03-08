a = 1234.57849; % your float point number
n = 16;         % number bits for integer part of your number      
m = 25;         % number bits for fraction part of your number
% binary number
d2b = fix(rem(a*pow2(-(n-1):m),2)); 
% the inverse transformation
b2d = d2b*pow2(n-1:-1:-m).'; 

addpath(genpath('.'));

clc;
p = Parameters();
pa = p.getParameters();

disp(pa.iTxAnt);
pa.iTxAnt = 6;

pu = p.getParameters();
disp(pu.iTxAnt);

keyboard;

lods = 4:8;%2:8;
sches = {'iterative','MRC','rotate','rotate2','noRelay','rotatesdr'};
%{'rotate','iterative'};%{'noRelay','MRC','novo'};
algs = 2;%4; [8 10]

for ss=1:length(sches)
    for ll=1:length(lods)
        %files = dir(['papero_*UE' num2str(lods(ll)) '*' sches{ss} '*.mat']);
        files = dir(['aftergre_*UE' num2str(lods(ll)) '*' sches{ss} '.mat']);
        
        %cdfs = figure;
        %hold all
        
        %for ff=1:length(files)
        disp(files.name)
        data = load(files.name);
        disp(data.sets{end})
        if iscell(data.ou.cellSNR)
            size(data.ou.cellSNR{1})
        else
            size(data.ou.cellSNR)
        end
    end
end

% adapt = false;
% 
% inter = 100;
% a = 0.0002;
% 
% users = 4;
% 
% %H = rand(4) + 1i*rand(4);
% 
% h = 2*(rand(4,1) + 1i*rand(4,1)) - 1 - 1i;
% h = h/norm(h);
% 
% hx = 2*(rand(4,1) + 1i*rand(4,1)) - 1 - 1i;
% hx = 2*hx/norm(hx);
% 
% w_ini = 2*(rand(4,1) + 1i*rand(4,1)) - 1 - 1i;
% w = w_ini/norm(w_ini);
% 
% %[abs(w'*h) abs(h'*h)]
% 
% %% -----------------------------------------------
% 
% Hu = [h hx]';
% iH = inv(Hu * Hu');
% 
% z3 = iH(1,2);
% s2 = pi - angle(z3);
% 
% w_t = (Hu')*iH*[1;exp(1i*s2)];
% 
% w_t = w_t/norm(w_t)
% 
% [abs(w_t'*h) abs(w_t'*hx)]
% 
% %% -------------------------
% 
% for iii=1:1
%     %% ---------
%     
%     d = w - h;
%     d = norm(d);
%     d = abs(d);
%     
%     a = (d^2)/2;
%     
%     % va = [a 0:0.0001:1];
%     vo = [];
%     %
%     % for aaa=1:length(va)
%     %
%     % a = va(aaa);
%     
%     % if d > sqrt(2)
%     %     disp('outro lado')
%     % else
%     h_ort = h - (w')*h*w;
%     
%     if norm(h_ort) == 0
%         disp('vixi!')
%         break;
%     end
%     
%     w_ort = h_ort/norm(h_ort);
%     
%     abs(w'*w_ort)
%     
%     if abs(w'*w_ort) > 0.005        
%         disp('vixi!')
%         break;
%     end
%     
%     dw = w_ort;
%     
%     
%     
%     %% projeta no plano w e w_ort
%     % 1 - a = w'*h
%     oa = 1 - abs(w'*h)
%     
%     
%     vu = linspace(0,oa,1000);
%     % vo = [];
%     %
%     for xu=1:length(vu)
%     %
%     % a = va(aaa);
%     
%     a = vu(xu);
%     
%     %a= oa;
%     
%     %beta = abs(dw'*h);
%     
%     beta = sqrt(2*a-a^2);
%     
%     ang = angle(w'*h) + angle(h'*w_ort);
%     
%     if angle(h'*w_ort) > 1e-8
%         keyboard;
%     end
%     
%     beta = beta*exp(1i*-ang);
%     
%     w_temp = (1-a)*w + beta*dw;
%     
%     w_temp = w_temp/norm(w_temp);
%     %end
%     
%     % vo(end+1) = abs(w_temp'*h);
%     %
%     % end
%     % [aux idx] = sort(vo,'descend');
%     % va(idx)
%     %w= w_temp;
%     
%     vo(end+1) = abs(w_temp'*h);
%     
%     end
%     
%     %[aux idx] = sort(vo,'descend');
%     %vu(idx)
%     
%     vo
%     vu
%     sum(diff(vo)<=0)
%     
%     w= w_temp;
%     
%     disp(['Corr: ' num2str(abs(w'*h))])
%     
% end

%% -------------------------------------------------

% snrs = zeros(1,inter);
% snrs2 = snrs;
% vtang = snrs;
% 
% snrs(1) = (abs(w'*h)^2); %/noisepw;
% vtang(1) = angle(w'*h);
% 
% conv = false;
% 
% for ii=2:inter
%     
%     h_ort = h - (w')*h*w;
%     
%     if norm(h_ort) == 0
%         break;
%     end
%     
%     %     if ii == 500
%     %         keyboard;
%     %     end
%     
%     %     vtang(ii) = norm(h_ort);
%     
%     w_ort = h_ort/norm(h_ort);
%     
%     if abs(w'*w_ort) > 0.5
%         break;
%     end
%     
%     dw = w_ort;
%     
%     beta = sqrt(2*a-a^2);
%     
%     ang = angle(w'*h) + angle(h'*w_ort);
%     
%     if angle(h'*w_ort) > 1e-8
%         keyboard;
%     end
%     
%     beta = beta*exp(1i*-ang);
%     
%     %     a1 = angle((beta*dw)'*h);
%     %     a2 = angle(((1-a)*w)'*h);
%     %
%     %     w = exp(1i*-(a1-a2))*w;
%     
%     w_temp = (1-a)*w + beta*dw;
%     
%     w_temp = w_temp/norm(w_temp);
%     
%     snrs(ii) = (abs(w_temp'*h)^2); %/noisepw;
%     
%     %     if snrs(ii) < snrs(ii-1)
%     %         ii
%     %     end
%     
%     if adapt
%         if snrs(ii) < snrs(ii-1)
%             if ii > 2
%                 if snrs(ii-2) < snrs(ii-1)
%                     a = a/2;
%                     conv = true;
%                 end
%             end
%         end
%     end
%     
%     if conv
%         %a2 = x;
%         snrs2(ii) = (abs(w_temp'*h)^2);
%     else
%         snrs2(ii) = (abs(w_temp'*h)^2);
%     end
%     
%     vtang(ii) = angle(w_temp'*h);
%     
%     w = w_temp;
% end
% 
% figure
% hold on
% plot(snrs,'xb')
% plot(vtang,'og')
% plot((abs(h'*h/norm(h))^2)*ones(1,inter),'sr')
% hold off
%
% H = rand(4) + 1i*rand(4);
%
% [vtM mtD] = matched(H,4,4);
%
%
% %[vtMo mtDo] = sdr(H,4,4);
%
%
% [vtMy mtDy] = umsf(H,4,4);
%
% [vtMn mtDn] = iurc(H,4,4);
%
% [vtMa mtDa] = iurcaux(H,4,4);
%
%
% w = rcc2(H);
% w = w/norm(w);
% w = isiu(H,w);
% %[min(abs(H*vtM)) min(abs(H*vtMo)) min(abs(H*vtMy))  min(abs(H*vtMn))]
% [min(abs(H*vtM)) min(abs(H*vtMy)) min(abs(H*vtMn)) min(abs(H*vtMa))]
% disp('end')