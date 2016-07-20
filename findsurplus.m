function [kspec,vspec,surplus] = findsurplus(grLs, grT, eval_no, tol)
% this function finds the first N=eval_no zeros of the secular function F(k)=det(1-exp(i*k*grLS).grT)
% these will be saved as kspec, and their corresponding eigenvectors as
% vspec
% the method used is finding all the eigenvalues of U(k)=exp(i*k*grLS).grT
% and propogate by k until one of them is 1, which gives a zero in the
% secular function.

    tic;
    
%defintions
    lmax = max(grLs);
    lmin = min(grLs);
    restol = tol*lmin;
    kspec = zeros(1,eval_no);
    surplus = zeros(1,eval_no);
    vspec = zeros(length(grT(1,:)),eval_no);
    
    count = 1;
    k = tol;
    
   %running the loop
    while (count < eval_no)
        %small check
        if k < kspec(count) + tol
            k = kspec(count)+tol;
            fprintf('k=%d got stuck in 1, k_num=%d\n',k,count);
        end
        U = diag(exp(1i*(k*grLs)))*grT;
        [V,D] = eig(U);
        evec = diag(D);
        %making the vector of angles 
        th_vec = angle(evec);
        %check if it is an eigenvalue and check for multiplicities
                d = sum(abs(th_vec) < restol);     
                if d==0
                        [yi,ii] = sort((th_vec+pi).*(th_vec<0));
                        ypos = min(th_vec(th_vec>0));               
                        delta_th = pi- yi(end-1);
                        Tmax = delta_th/lmax;
                        k_next = k + Tmax;

                        if pi- yi(end) < lmin*Tmax
                            %we have an eigenvalue in the interval...go
                            %Newton!
                                        eigenvec = V(:,ii(end));
                                        tprime = real((eigenvec'.*grLs)*eigenvec);
                                        dk = (pi-yi(end))/tprime;
                                        dk = min(dk, max(0,(ypos + pi - yi(end))/(lmax-lmin)));
                                     %Newton
                                        k_in = k;

                                        [k_out,v] = Newton(grLs, grT, k_in, dk,tol);
                                                %cheking for multiplicity
                                                if k_out > k_next
                                                    k = k_out;
                                                    fprintf('k=%d may have multiplicity, k_num=%d \n',k,count);
                                                    continue
                                                end
                                                %checking that we are not going
                                                %backwards
                                                if k_out < kspec(count) + tol
                                                     k = max(kspec(count),k_in)+2*tol;
                                                     fprintf('k=%d got stuck in 2,k_num=%d \n',k,count);
                                                     continue
                                                end
                                          
                                        count = count+1;
                                        kspec(count) = k_out;
                                        vspec(:,count) = v;
                                        surplus(count) = nodalcount(grLs,vspec(:,count),kspec(count))-(count-1);
                                        k =max(k_next,k_out+tol);
                                        continue
                                                                                         
                         else
                            k = k + (pi-yi(end))/lmax;
                            continue
                        end
                %if we got multiplicity!        
                else
                    fprintf('k=%d is of degree %d,k_num=%d \n',k,d,count);

                    [~,id] = sort(pi-abs(th_vec));
                    for j=1:d
                        count = count+1;
                        kspec(count) = k;
                        vspec(:,count) = V(:,id(end-j+1));
                        surplus(count) = nodalcount(grLs,vspec(:,count),kspec(count))-(count-1);
                    end
                    k = kspec(count) + tol;
                end
    end
        
       
        
    
        
        toc;
