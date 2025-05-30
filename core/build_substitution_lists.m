function [subs_vars, subs_vals] = build_substitution_lists(x, dN, x_num, dN_num, Dx, Dx_num)
%BUILD_SUBSTITUTION_LISTS Create substitution pairs for symbolic x, dN, and Dx
%   [subs_vars, subs_vals] = build_substitution_lists(x, dN, x_num, dN_num, Dx, Dx_num)
%   returns cell arrays for substituting numerical values into symbolic x, dN, and Dx

  % Get sizes
  [nk, nd]      = size(x);
  [~, nxi, nj]  = size(dN);

  subs_vars = {};
  subs_vals = {};

  % x substitution
  for k = 1:nk
    for d = 1:nd
      subs_vars{end+1} = x(k,d);
      subs_vals{end+1} = x_num(k,d);
    end
  end

  % dN substitution
  for k = 1:nk
    for xi = 1:nxi
      for j = 1:nj
        subs_vars{end+1} = dN(k,xi,j);
        subs_vals{end+1} = dN_num(k,xi,j);
      end
    end
  end

  % Dx substitution (if provided)
  if nargin >= 6 && ~isempty(Dx) && ~isempty(Dx_num)
    for k = 1:nk
      for d = 1:nd
        subs_vars{end+1} = Dx(k,d);
        subs_vals{end+1} = Dx_num(k,d);
      end
    end
  end

end
