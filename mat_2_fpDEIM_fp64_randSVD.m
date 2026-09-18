
    f1 = figure('Visible','off');
    ax1 = axes(f1);

    f2 = figure('Visible','off');
    ax2 = axes(f2);

    err = [];
    err_eta_q = [];
    err_eta_p = [];
    
    [~,S,~]    = svd(full(A),0); S = diag(S);
    [Vs,~,Ws]  = rsvd_double(A,maxK);
 
    for kk = 1:length(formats)
        options.format = formats(kk);
        err = [];
        err_eta_q = [];
        err_eta_p = [];

        for k=1:maxK 
            
            [C,U,R,eta_p,eta_q] = CUR(A,Vs(:,1:k),Ws(:,1:k),options);
            err = [err; norm(full(A-C*U*R),2)];
            err_eta_q = [err_eta_q; eta_q];
            err_eta_p = [err_eta_p; eta_p];

        end
        semilogy(ax1,1:size(err,1),err, 'Marker', markerShapes{kk},'Color',colorShapes{kk});
        hold(ax1, 'on');
        legendLabels{kk} = sprintf("$u_D$ = %s",formats_full(kk));
        semilogy(ax2,1:size(err_eta_p,1),err_eta_p, 'Marker', markerShapes{kk},'Color',colorShapes{kk});
        hold(ax2, 'on');
    end


    semilogy(ax1,1:maxK,S(2:maxK+1),'ks-');
    hold(ax1, 'on');
    legendLabels{end+1} = '$\sigma_{k+1}$';

    hold(ax1, 'off');
    grid(ax1,'on');
    legend(ax1,legendLabels,'Location', 'best','Interpreter', 'latex','FontSize',14)
    xlabel(ax1,'k','FontSize',18);
    ylabel(ax1,'$$||A-\hat{C} U \hat{R}||$$', 'Interpreter', 'latex','FontSize',18);
    
        title(ax1,'Finite precision DEIM for fp64 randSVD','Interpreter', 'latex','FontSize',18);

    hold(ax2, 'off');
    grid(ax2,'on');
    legend(ax2,legendLabels,'Location', 'best','Interpreter', 'latex','FontSize',14)
    xlabel(ax2,'k','FontSize',18);
    ylabel(ax2,'$$\hat{\eta}_p$$', 'Interpreter', 'latex','FontSize',18);
    
        title(ax2,'$\hat{\eta}_p$ for finite precision DEIM for fp64 randSVD', 'Interpreter', 'latex','FontSize',18);

saveas(f1,'mat_2_fpDEIM_fp64_randSVD','pdf')
saveas(f2,'mat_2_eta_DEIM_fp64_randSVD','pdf')

return 

