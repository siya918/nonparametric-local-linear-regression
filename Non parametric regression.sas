options ls=72 nocenter;
dm 'odsresults;clear';
title;

/* Generate dataset */
data work.npdat2017;
    call streaminit(42);
    do i = 1 to 1000;
        x = rand('uniform') * 10;
        y = 3 * sin(x) + 0.5 * x + rand('normal', 0, 1.2);
        output;
    end;
    drop i;
run;

proc sort data=work.npdat2017;
    by x;
run;

/* Graph 1 - Raw scatter plot */
proc gplot data=work.npdat2017;
    plot y*x;
run;
quit;

/* Local linear regression */
proc iml;
    use work.npdat2017;
    read all into xy;
    x1 = xy[,1];
    y1 = xy[,2];
    n  = nrow(xy);
    m  = 30;
    print n m;

    yh = j(0,1,.);

    do i = 1 to n;
        f_point = x1[i,];
        x0      = x1 - f_point;
        axy     = abs(x0) || x0 || xy;
        call sort(axy, {1});
        x       = j(m,1,1) || axy[1:m,2];
        y       = axy[1:m,4];
        bh      = inv(x`*x) * x`*y;
        yh_fp   = bh[1,1];
        yh      = yh // yh_fp;
    end;

    res = xy || yh;
    nm  = {"x" "y" "yh"};
    print res[colname=nm];

    create work.npres from res[colname=nm];
    append from res;
quit;

/* Graph 2 - Observed vs fitted overlay */
symbol1 interpol=none width=4
        color=blue
        value=dot
        height=2;

symbol2 interpol=line width=3
        color=red
        value=dot
        height=.3;

proc gplot data=work.npres;
    plot (y yh)*x / overlay;
run;
quit;

proc iml;
    use work.npdat2017;
    read all into xy;
    x1 = xy[,1];
    y1 = xy[,2];
    n  = nrow(xy);

    /* candidate bandwidth values */
    m_cand = do(5, 55, 5);
    n_cand = ncol(m_cand);
    cv_res = j(n_cand, 2, .);

    do k = 1 to n_cand;
        m        = m_cand[k];
        loo_mse  = 0;

        do i = 1 to n;
            /* remove observation i */
            idx_train = (1:n)[loc((1:n) ^= i)];
            x_train   = x1[idx_train];
            y_train   = y1[idx_train];
            x_test    = x1[i];
            y_test    = y1[i];

            /* distances from test point */
            x0  = x_train - x_test;
            axy = abs(x0) || x0 || x_train || y_train;
            call sort(axy, {1});

            /* m nearest neighbours */
            x_nn = j(m,1,1) || axy[1:m, 2];
            y_nn = axy[1:m, 4];

            /* local OLS */
            bh   = inv(x_nn`*x_nn) * x_nn`*y_nn;
            yh_i = bh[1,1];

            loo_mse = loo_mse + (y_test - yh_i)**2;
        end;

        cv_res[k,1] = m;
        cv_res[k,2] = loo_mse / n;
    end;

    nm = {"m" "LOOCV_MSE"};
    print cv_res[colname=nm];

    /* find optimal m */
    min_idx   = cv_res[,2][<:>];
    optimal_m = cv_res[min_idx, 1];
    print optimal_m;

    create work.cvres from cv_res[colname=nm];
    append from cv_res;
quit;

/* Plot LOOCV error vs bandwidth */
symbol1 interpol=line width=3
        color=steelblue
        value=dot
        height=1.5;

proc gplot data=work.cvres;
    plot loocv_mse*m;
    title 'LOOCV MSE vs Bandwidth';
run;
quit;