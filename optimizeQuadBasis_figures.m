function optimizeQuadBasis_figures(quad,threshold,maxSwaps)
% given a quadBasis representing X in factored form, reduces it so that all
% elements of Xopt are below a threshold by working on factors and produces
% the figures for getFiguresForOptimizingRandomQuad
%
% A version of optimizeQuadBasis that returns the sequence of |det(X)| and
% max |X| along the steps, for plotting.
%
% (c) 2015-2015 F. Poloni <poloni@math.tu-berlin.de> and others
% see AUTHORS.txt and COPYING.txt for details
% https://bitbucket.org/fph/pgdoubling

if not(exist('threshold','var')) || isempty(threshold)
    threshold=2;
end

if threshold<1+sqrt(eps(class(quad.X)))
    error('PGDoubling:thresholdTooSmall','you can only hope to enforce thresholds S=1+sqrt(eps) or larger');
end

if not(exist('maxSwaps','var')) || isempty(maxSwaps)
    maxSwaps=20*size(quad.X,2);
end

absdets = [];
maxnorms = [];
%optimization
swaps=0;
invcond=1;
indices = 1:length(quad.X); %helper for Matlab indexing
iterations = 0;

figure
while(true)
    iterations = iterations + 1;
    sym = symBasisFromQuadBasis(quad);
    absdets = [absdets abs(det(sym.X))];
    maxnorms = [maxnorms max(max(abs(sym.X)))];
    
    switch iterations
        case 1
            sp(1) = subplot(2,2,1);imagesc(abs(sym.X),[0 2*threshold]);
            title('(a) iteration 1','fontweight','normal','fontsize', 11);
        case 10
            sp(2) = subplot(2,2,2);imagesc(abs(sym.X),[0 2*threshold]);
            title('(b) iteration 10','fontweight','normal','fontsize', 11);
        case 20
            sp(3) = subplot(2,2,3);imagesc(abs(sym.X),[0 2*threshold]);
            title('(c) iteration 20','fontweight','normal','fontsize', 11);
    end
    % TODO: use Natasa's heuristic instead of recomputing norms each time
    % column norms of C
    squaredColNorms = sum(abs(quad.X(quad.v,quad.v)).^2,1);
    [maxSquaredColNorm,maxi] = max(squaredColNorms);
    if maxSquaredColNorm > threshold
        % I can't find how to avoid this index juggling --federico
        indicesInV = indices(quad.v);
        maxi = indicesInV(maxi);
        [quad.X, quad.v, stepcond] = updateQuadBasisOut(quad.X, quad.v, maxi);
        swaps = swaps + 1;
    else
        % row norms of B
        squaredRowNorms = sum(abs(quad.X(~quad.v,~quad.v)).^2,2);
        [maxSquaredRowNorm,maxj] = max(squaredRowNorms);
        if maxSquaredRowNorm > threshold
            indicesNotInV = indices(~quad.v);
            maxj = indicesNotInV(maxj);
            [quad.X, quad.v, stepcond] = updateQuadBasisIn(quad.X, quad.v, maxj);
            swaps = swaps + 1;
        else
            [maxvec, maxis] = max(abs(quad.X(~quad.v,quad.v)));
            [maxval, maxj] = max(maxvec);
            if maxval > threshold
                indicesNotInV = indices(~quad.v);
                maxi = indicesNotInV(maxis(maxj));
                indicesInV = indices(quad.v);
                maxj = indicesInV(maxj);
                [quad.X, quad.v, stepcond] = updateQuadBasisInOut(quad.X, quad.v, maxi, maxj);
                swaps = swaps + 2;
            else
                break %finished optimization
            end
        end
    end
    invcond = invcond*stepcond;
    if swaps >= maxSwaps
        warning('PGDoubling:stagnated','failed to produce a X with elements below the required threshold. Try running with a larger threshold.');
        break
    end
end

% add the final snapshot
  sp(4) = subplot(2,2,4);imagesc(abs(sym.X),[0 2*threshold]);
  title(['(d) iteration ',int2str(iterations)],'fontweight','normal','fontsize', 11);
  
    
%add the joint colorbar on the right of plots
h = colorbar;
set(h, 'Position', [.9 .11 .05 .8150]);
for i=1:4
pos=get(sp(i), 'Position');
set(sp(i), 'Position', [pos(1) pos(2) 0.8*pos(3) pos(4)]);
end;

%set(gcf,'Color','w');
%export_fig ../figs/snapshots.pdf
%export_fig ../figs/snapshots.eps

figure
    options = {'Interpreter','latex'};
    subplot(1,2,1);plot(maxnorms,'k');xlabel('iteration','fontsize',14);ylabel('$\max |X_{ij}|$', options{:},'fontsize',14,'rotation',0,'fontsize',14,'verticalalignment','top','horizontalalignment','right');
    %axis([1 iterations 0 max(maxnorms)]);
    axis image;
    subplot(1,2,2);semilogy(absdets,'k');xlabel('iteration','fontsize',14);ylabel('$|\det X|$', options{:},'fontsize',14,'rotation',0,'fontsize',14,'verticalalignment','top','horizontalalignment','right');
    %axis([1 iterations 0 max(absdets)]);
    axis tight;

    %set(gcf,'Color','w');
    %export_fig ../figs/iterations.pdf
    %export_fig ../figs/iterations.eps

    % rescale the figure produced before exporting
    
sum(quad.v')
fprintf(1,'max(|XOpt_{ij}|) = %f\n',maxnorms(iterations))



