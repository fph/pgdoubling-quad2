function oldQuad = quadBasisFromQuadBasis(oldQuad,newv)
% Transform oldQuad.X into a newv-semidefinite matrix by working directly
% with the blocks; section 5 from the paper.

oldv = oldQuad.v; % oldv defines \mathcal{I}, newv defines \mathcal{J}
vPPT = logical(xor(oldv,newv)); % defines \mathcal{K} for the PPT

% Determine which case we are in for the PPT
if all(ismember(find(newv),find(oldv))) % true iff \mathcal{J} \subseteq \mathcal{I}; Case 1
    vPPTlocal = vPPT(oldv);
    [~,oldQuad.X(oldv,oldv)] = leftUnitary(oldQuad.X(oldv,oldv),~vPPTlocal,vPPTlocal);
    inverse = inv(oldQuad.X(vPPT,vPPT));
    oldQuad.X(vPPT,vPPT) = inverse;
    oldQuad.X(vPPT,newv) = -inverse*oldQuad.X(vPPT,newv);
    oldQuad.X(~oldv,newv) = oldQuad.X(~oldv,newv) + oldQuad.X(~oldv,vPPT)*oldQuad.X(vPPT,newv);
    oldQuad.X(~oldv,vPPT) = oldQuad.X(~oldv,vPPT)*inverse;
    oldQuad.X(vPPT,~oldv) = 0;
else if all(ismember(find(oldv),find(newv))) % true iff \mathcal{J} \supseteq \mathcal{I}; Case 2
            vPPTlocal = vPPT(~oldv);
            [~,oldQuad.X(~oldv,~oldv)] = rightUnitary(oldQuad.X(~oldv,~oldv),vPPTlocal,~vPPTlocal);
            inverse = inv(oldQuad.X(vPPT,vPPT));
            oldQuad.X(oldv,vPPT) = 0;
            oldQuad.X(vPPT,vPPT) = inverse;
            oldQuad.X(~newv,vPPT) = oldQuad.X(~newv,vPPT)*inverse;
            oldQuad.X(~newv,oldv) = oldQuad.X(~newv,oldv) - oldQuad.X(~newv,vPPT)*oldQuad.X(vPPT,oldv);
            oldQuad.X(vPPT,oldv) = -inverse*oldQuad.X(vPPT,oldv);
    else % Case 3
        % Blocks for the pivot matrix globally
        v = logical(oldv.*vPPT); % \mathcal{I} \cap \mathcal{K}
        w = logical(~oldv.*vPPT); % \mathcal{I}^c \cap \mathcal{K}
        u = logical(oldv.*(~vPPT)); % \mathcal{I} \cap \mathcal{K}^c
        z = logical(~oldv.*(~vPPT)); % \mathcal{I}^c \cap \mathcal{K}^c
            
        % Blocks for leftUnitary and rightUnitary locally
        vlocal = v(oldv);
        wlocal = w(~oldv);
        [~,oldQuad.X(oldv,oldv)] = leftUnitary(oldQuad.X(oldv,oldv),~vlocal,vlocal);
        [~,oldQuad.X(~oldv,~oldv)] = rightUnitary(oldQuad.X(~oldv,~oldv),wlocal,~wlocal);
        
        [N,K,L,H] = inverseQuasidef(oldQuad.X(v,v),oldQuad.X(w,v),oldQuad.X(w,w));

        m = length(L);
        C21new = H(m+1:end,1:m)'*oldQuad.X(v,u)-L*oldQuad.X(w,u);
        B21new = oldQuad.X(z,w)*H(1:m,m+1:end)+oldQuad.X(z,v)*N;
        A11new = -N*H(m+1:end,m+1:end)'*oldQuad.X(v,u) - K*oldQuad.X(w,u);
        A21new = oldQuad.X(z,u) + oldQuad.X(z,v)*A11new + oldQuad.X(z,w)*oldQuad.X(w,w)'*K'*oldQuad.X(v,v)'*oldQuad.X(v,u) - oldQuad.X(z,w)*H(1:m,1:m)*L*oldQuad.X(w,u);
        A22new = oldQuad.X(z,v)*K+oldQuad.X(z,w)*H(1:m,1:m)*L;
        
        oldQuad.X(u,w) = 0;
        oldQuad.X(w,u) = C21new;
        oldQuad.X(w,w) = L;
        oldQuad.X(v,u) = A11new;
        oldQuad.X(v,w) = K;
        oldQuad.X(v,v) = N;
        oldQuad.X(v,z) = 0;
        oldQuad.X(z,u) = A21new;
        oldQuad.X(z,w) = A22new;
        oldQuad.X(z,v) = B21new;
    end
end

oldQuad.v=newv;

end

