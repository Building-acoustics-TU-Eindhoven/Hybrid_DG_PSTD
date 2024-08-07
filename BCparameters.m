classdef BCparameters
% class containing a pole residue model of the reflection coefficient

    properties

        Kind   % 0==no poles, 1==only real poles, 2==complex poles
        BCcode % boundary code associated with pole residue model
        Rinf   % direct term
        C      % real residues
        A      % real part of complex conjugate residues
        B      % imag part of complex conjugate residues
        lambda % real poles (lambda < 0)
        alpha  % real part of complex conjugate poles (alpha < 0)
        beta   % imag part of complex conjugate poles
               % --------------------------------------------------------------
               % COMPLEX CONJUGATE PAIRS CORRESPOND TO ONE ENTRY (i.e., the
               % complex conjugate of each entry is not explicitly defined)
        R      % frequency response of the reflection coefficient
        poles
        residues

    end %properties

    methods

        function obj = BCparameters(BCcode,Rinf,poles,residues)
            if nargin==2
                poles = [];
                residues = [];
            end

            assert(all(size(poles)==size(residues)), ...
                   'different number of poles and residues')

            obj.BCcode = BCcode;
            obj.Rinf = Rinf;
            obj.Kind = 0;
            obj.poles = poles;
            obj.residues = residues;

            obj.R = @(s) obj.Rinf;

            % if poles/residues are defined
            if length(poles)
                % if no imaginary part
                if all(imag(poles)<eps & imag(residues)<eps)
                    % only real poles
                    obj.Kind = 1;
                    obj.C = real(residues);
                    obj.lambda = real(poles);
                    obj.R = @(s) obj.R(s) + sum(obj.C./(s-obj.lambda));
                else % all the poles are considered complex
                    obj.Kind = 2;
                    % fix amplitude of real poles/residues (since they
                    % are considered complex)
                    idx = (imag(poles)<eps & imag(residues)<eps);
                    residues(idx) = residues(idx)/2;
                    obj.A = real(residues);
                    obj.B = imag(residues);
                    obj.alpha = real(poles);
                    obj.beta = abs(imag(poles)); % make imag part positive, so
                                                 % that they all have the same sign
                    obj.R = @(s) obj.R(s) + sum( ...
                        (obj.A+1i*obj.B)./(s-(obj.alpha+1i*obj.beta)) + ...
                        (obj.A-1i*obj.B)./(s-(obj.alpha-1i*obj.beta)) ...
                        );
                end
            end

        end



    end %methods

    methods (Static)

        function [BCcodesInst,BCcodesReal,BCcodesCplx]=sort(BCcodes,BCparams)
        % sort the pole-residue models according to their kind

            BCcodesInst = [];
            BCcodesReal = [];
            BCcodesCplx = [];

            for ii=1:length(BCparams)

                % find position of corresponding BCcode
                ibc = find(BCcodes==BCparams{ii}.BCcode);

                % safety checks
                if ~length(ibc), error('boundary code not found!'); end
                if length(ibc)>1, error('boundary code defined multiple times!'); end

                switch BCparams{ii}.Kind
                  case 0
                    BCcodesInst = [BCcodesInst ibc];
                  case 1
                    BCcodesReal = [BCcodesReal ibc];
                  case 2
                    BCcodesCplx = [BCcodesCplx ibc];
                end

            end
        end %sort

    end %static methods

end %class

% $$$ function [poles_real, resis_real, poles_cplx, resis_cplx] = process_pr(poles,residues)
% $$$ % function to separate the poles and residues into real p/r and complex p/r
% $$$
% $$$     poles_real = [];
% $$$     resis_real = [];
% $$$     poles_cplx = [];
% $$$     resis_cplx = [];
% $$$
% $$$     for ii=1:length(poles)
% $$$         % if no imaginary part
% $$$         if imag(poles(ii))<eps && imag(residues(ii))<eps
% $$$             % assume real pole and residue
% $$$             poles_real = [poles_real real(poles(ii))];
% $$$             resis_real = [resis_real real(residues(ii))];
% $$$         else % complex pole and residue
% $$$             poles_cplx = [poles_cplx poles(ii)];
% $$$             resis_cplx = [resis_cplx residues(ii)];
% $$$         end
% $$$     end
% $$$
% $$$ end
