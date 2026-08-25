function R = HB_residual_TLCD_2GDL(X,system,H,N)
% Extrai parâmetros
M = system.M;
D = system.D;
K = system.K;
Fex1 = system.Fex1;
epsilon = system.epsilon;
L = system.L;
n = system.n;

% Conversion of sine-cosine to complex-exponential representation

I0 = 1:n; % Indice dos Coef. constantes para cada GDL
ID = n+(1:H*n); % Indice correspondente aos harmonicos de 1 até H  

IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); % Indice indicando os Coef. de cossenos
IS = IC+n; % Indice dos Coef. de senos

Q = zeros(n*(H+1),1); % Montando a estrutura dos Coef. de Fourier Q
Q(I0) = X(I0); % Copiando a componente do coef. constante de X para Q
Q(ID) = X(IC)-1i*X(IS); %Conversão cos/sen para a representação complexa

% Setup excitation vector
Fex = zeros(n*(H+1),1); % Montando a estruturo do vetor excitação
h = 1;% Forcing of h-th harmonic (O forçamento acontece na componente do cos(Omega*t))
Fex(h*n+(1:n)) = Fex1; % Adicionando a excitação

% Excitation frequency (Frequência de excitação)
Om = X(end);

[Fnl] = HB_nonlinear_forces_AFT(Q,Om,H,N,...
    system.nonlinear_elements);

% Dynamic force equilibrium

Rc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*Q + Fnl - Fex;

% Conversion from complex-exponential to sine-cosine representation
R = zeros(size(X,1)-1,1);
R(I0) = real(Rc(I0));
R(IC) = real(Rc(ID));
R(IS) = -imag(Rc(ID));
%%

function [F] = ...
    HB_nonlinear_forces_AFT(Q,Om,H,N,nonlinear_elements)

%Criando o vetor das forças não-lineares 
F = zeros(size(Q)); 

%% Iteração dos elementos não-lineares
for nl=1:length(nonlinear_elements)
    % Especificando amostras de tempo ao longo do período
    tau = (0:2*pi/N:2*pi-2*pi/N)';
    
    % Determinando a direção da força não-linear 
    w = nonlinear_elements{nl}.force_direction;
    W = kron(eye(H+1),w);
    
    % Convertendo para o Domínio do Tempo 
    %(Aplicando o inverse discrete Fourier transform)
    H_iDFT = exp(1i*tau*(0:H));
    % Solução adotada
    qnl = real(H_iDFT*(W'*Q)); % Deslocamento
    qnldot = real(H_iDFT*(1i*Om*(0:H)'.*(W'*Q))); % Velocidade
    
    % Avaliando a Força NL no Domínio do Tempo
    switch lower(nonlinear_elements{nl}.type)
        case 'cubicspring'
            fnl = nonlinear_elements{nl}.stiffness*qnl.^3;
        
        case 'quadraticdamper'
            fnl = nonlinear_elements{nl}.damping*abs(qnldot).*qnldot;
    end
    %% Calculada as forças no Domínio do Tempo, precisamos voltar para o 
    %domínio da frequência
    
    % Aplicando a FFT
    Fnlc = fft(fnl(end-N+1:end))/N;
    
    
    % A linha extrai o termo constante (índice 1)
    % e os próximos H harmônicos (índice 2 até H+1)
    Fnl = [real(Fnlc(1));2*Fnlc(2:H+1)];
    
    F = F + W*Fnl;
end
