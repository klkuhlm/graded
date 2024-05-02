
bores = {'DBT10','DBT11','DBT12','DBT13','L4B01'};
for j = [1,2,3,4,5]

file= bores{j};

%load('DREAM_ZS.mat');
type = 'powerlaw'; 
load(['powerlaw_',file,'_results.mat']);

nchains = size(Sequences,3);
npar = size(Sequences,2)-2;
niter = size(Sequences,1); % iterations per chain
%%niter = 762000/nchains;

best_fit = zeros(npar,1);

% pflotran vars
vars = {'k_0','n_0','\eta','\tau','c_m','m'}; % ,'\eta','\tau'

colors = {'r','g','b','k','c','m','r','g','b','c','m','k'};
title = [file,' brine flux'];

% plot of chains
% ***************
if 1
    f = figure();
    %f.Renderer = 'Painters';
    ha = tight_subplot(npar,1,.005,[.1 .01],[.1 .01]);
    for i = 1:npar
        %subplot(npar,1,i);
        % plot each chain as a differnt colors
        for j =1:nchains
            plot(ha(i),Sequences(:,i,j),[colors{j} '.'],'MarkerSize',0.5);
            hold(ha(i),'on');
        end
        ylabel(ha(i),vars{i},'FontSize',12);
        xlim(ha(i),[0,niter]);
        if i < npar
            set(ha(i),'XTickLabel',[]);
        end
    end

    xlabel(ha(npar),'iterations per chain','FontSize',9);
    suptitle(title);
    f.PaperUnits = 'inches';
    f.PaperSize = [6 5];
    f.PaperUnits = 'normalized';
    f.PaperPosition = [0 0 1 1];
    print([file,'-all-chains.png'],'-dpng','-r200');
    %savefig(f,[file,'-all-chains.fig']);
end

% reshape output after burnin period
burnin = 50000;
jump = niter-burnin;
tmp = zeros(npar+1,jump*nchains,1);
for i = 1:npar
    for j=1:nchains
       tmp(i,jump*(j-1)+1:j*jump) = Sequences(burnin+1:niter,i,j);
    end
end

tmp(tmp == 0.0) = NaN;

% plot of histograms
% ******************

nbins = 25;

plothist = 1;

if plothist
  f = figure();
  %f.Renderer = 'Painters';
end

for i = 1:npar
    subplot(2,3,i);
    % normalized histogram  from 
    % http://www.mathworks.com/matlabcentral/fileexchange/
    % 22802-normalized-histogram
    xx = tmp(i,:);
    [n,xout] = histnorm(xx,nbins);
    m = mean(xx(~isnan(xx)));
    v = var(xx(~isnan(xx)));
    if plothist
        bar(xout,n,'BarWidth',1.0);
        hold on;
        xlabel(vars{i},'FontSize',12);
        ylabel('Posterior density','FontSize',9);

        % normal dist
        plot(xout,exp(-(xout-m).^2 ./(2*v))./sqrt(2*pi*v),'r-');
        grid on;
    end
    [mx,ix] = max(n); % estimate mode
    fprintf("%s\tmean=%d\tvar=%d\tmode=%d\n",vars{i},m,v,xout(ix));
    best_fit(i) = xout(ix); % mode of distribution
end


if plothist
    suptitle(title);
    f.PaperUnits = 'inches';
    f.PaperSize = [7,5];
    f.PaperUnits = 'normalized';
    f.PaperPosition = [0 0 1 1];
    print([file,'-posterior-histograms.png'],'-dpng','-r200');
    %savefig(f,[file,'-posterior-histograms.fig']);
end


% plot of correlation between posterior distributions (scatter)
% ***************************************************
if 0
  f = figure();
  %f.Renderer = 'Painters';
  [H,AX,BigAx,P,PAx] = plotmatrix(transpose(tmp(1:npar,:)),'r.');
  for i = 1:npar
     ylabel(AX(i,1),vars{i},'FontSize',9);
     xlabel(AX(npar,i),vars{i},'FontSize',9);
  end
  suptitle(title);
  f.PaperUnits = 'inches';
  f.PaperSize = [6,6];
  f.PaperUnits = 'normalized';
  f.PaperPosition = [0 0 1 1];
  print([file,'-posterior-joint-distributions.png'],'-dpng');
end

% plot of correlation between posterior distributions (2d colors)
% ***************************************************
if 1
  f = figure();
  %f.Renderer = 'Painters';
  f.PaperUnits = 'inches';
  f.PaperSize = [6,6];
  %suptitle(sample_names{k});
  ha = tight_subplot(npar,npar,.005,[.1 .01],[.1 .01]);
  for i = 1:npar % rows
      for j = 1:npar % cols
          idx = npar*(i-1)+(j-1)+1;
          if i > j || j > i
              histogram2(ha(idx),tmp(j,:),tmp(i,:),...
                  'DisplayStyle','tile','EdgeColor','none');
              if i < npar
                  set(ha(idx),'Xtick',[]);
              end
              if j > 1
                  set(ha(idx),'Ytick',[]);
              end
          elseif i == j
              xx = tmp(i,:);
              [n,xout] = histnorm(xx,nbins);
              bar(ha(idx),xout,n,'BarWidth',1.0);
              %hold on
              %m = mean(xx(~isnan(xx)));
              %v = var(xx(~isnan(xx)));
              % normal dist
              %plot(ha(idx),xout,exp(-(xout-m).^2 ./(2*v))./sqrt(2*pi*v),'r-');
              grid on;
          end
          if i < npar
              set(ha(idx),'Xtick',[]);
          else
              xlabel(ha(idx),vars{j},'FontSize',9);    
          end
          if j > 1
              set(ha(idx),'Ytick',[]);
          else             
              ylabel(ha(idx),vars{i},'FontSize',9);
          end
          grid(ha(idx),'on');
      end
  end
  
  f.PaperUnits = 'normalized';
  f.PaperPosition = [0 0 1 1];
  tightfig; 
  print([file,'-posterior-joint-histogram2d.png'],'-dpng','-r200');
  %savefig(f,[file,'-posterior-joint-histogram2d.fig']);
end  


% horsetail plot (must be re-computed)
% *******************************************************
if 1
    nhorsetail = 100; % samples per chain
    f = figure();
    %f.Renderer = 'Painters';  % needed to make vector figure
           
    for i = 1:nhorsetail
        for j = 1:nchains
           
            pars = Sequences(niter-i,1:npar,j);
            evalstr = ['Qmod = ',type,'(pars,Extra);' ];
            eval(evalstr);
                                                                
            subplot(1,2,1); % linear
            p = plot(Extra.time,Qmod(:),'r-');
            p.Color = [0.5,0.5,0.5];
            hold('on');
            subplot(1,2,2); % log
            p = loglog(Extra.time,Qmod(:),'r-');
            p.Color = [0.5,0.5,0.5];
            hold('on');
        end
    end

    subplot(1,2,2);
    loglog(Extra.time,Measurement.MeasData(:),'ro');
    subplot(1,2,1);
    plot(Extra.time,Measurement.MeasData(:),'ro');
    
    suptitle(title);
    f.PaperUnits = 'inches';
    f.PaperSize = [6,6];
    f.PaperUnits = 'normalized';
    f.PaperPosition = [0 0 1 1];
    print([file,'-horsetail.png'],'-dpng','-r300');
    %print([file,'-horsetail.eps'],'-depsc');
end
close all;
end
