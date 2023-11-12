using Plots, LinearAlgebra,Distributions

#Recomendaciones generales para mejorar el CSD-3
#Fuente de luz extendida
#f(x)=x^2 
#f(x,y)=x^2+y^2
#Se puede nombrar dos operaciones con el mismo nombre de una función 
#f.(x,y) Broadcast -> operar vectores
#broadcast (f.())
#f.(x,y') y' adjunta para crear una matriz de arrays 

bits=512; #Cantidad de puntos del "Meshgrid"
N=20; #Número de puntos random/Número de cortes transversales/Número de vórtices 

w0=.5; # Cintura del haz [mm]
kx=10; #Vector de propagación en (x,y) (Da la inclinación del vector)
a=.5*w0; #Tamaño de la fuente
A=1;
l=1;

function U(X,Y,w0,kx)
    u= [exp.(-(x.^2+y.^2)./(w0^2)) + exp.(kx*im*y).*exp.(-(x.^2+y.^2)./(w0^2)) for x in X, y in Y] 
   return abs.(u)
end

function Vortex(X,Y,w0,A,l)
    u=[A/(w0^l).*(sqrt.(x.^2+y.^2).^abs(l)).*exp.(-(x.^2+y.^2)./(w0^2)).*exp.(im.*l.*atan.(y,x)) for x in X,y in Y]
    return u
end

function I_Gauss(w0,kx,a,bits,N)

    xmax=2*w0; # Tamaño de la rendija

    UU=zeros(bits,bits,N);

    X = range(-xmax,xmax*(bits-2)/bits,length=bits); 
    Y = range(-xmax,xmax*(bits-2)/bits,length=bits);

    r=a*sqrt.(rand(1,N)); #Puntos random en la coordenada radial
    θ=2*pi*rand(1,N); #Puntos random en la coordenada radial

    rx=r.*cos.(θ); 
    ry=r.*sin.(θ);

    for i=1:N
        UU[:,:,i]=U(X .- rx[i],Y .- ry[i],w0,kx) 
    end

    I=mean(UU,dims=3);

    heatmap(X,Y,I[:,:],colormap=:grays)

end

function I_X_Vortex(w0,a,A,l,bits,N)

    xmax=2*w0; # Tamaño de la rendija

    UU=zeros(Complex,bits,bits,N); #Complex -> Tipo de valor de los argumentos de la matriz 

    X = range(-xmax,xmax*(bits-2)/bits,length=bits); #Collect= vector 
    Y = range(-xmax,xmax*(bits-2)/bits,length=bits);

    r=a*sqrt.(rand(1,N)); #Puntos random en la coordenada radial
    θ=2*pi*rand(1,N); #Puntos random en la coordenada radial

    rx=r.*cos.(θ) 
    ry=r.*sin.(θ) 

    for i=1:N
        UU[:,:,i]=Vortex(X .- rx[i],Y .- ry[i],w0,A,l) 
    end

    (abs.(UU)).^2

    In=mean((abs.(UU)).^2,dims=3); #Intensidad 

    Xn=mean(real.(UU.*(-UU)),dims=3) #Correlación cruzada

    a_1=heatmap(X,Y,In[:,:],colormap=:grays,aspect_ratio=1,title="⟨I(r)⟩",titlefont=font(12,"Computer Modern")) #Aspect ratio -> Centrar la visualización 
    a_2=heatmap(X,Y,Xn[:,:],colormap=:grays,aspect_ratio=1,title="⟨X(r)⟩",titlefont=font(12,"Computer Modern"))

    plot(a_1,a_2,plot_title="Vórtice óptico [N=$N A=$A w0=$w0 a=$a l=$l]",plot_titlevspan=.09,plot_titlefont=font(14,"Computer Modern"))

end

#a variable, l fijo 

a=2*w0;
A_2=I_X_Vortex(w0,a,A,l,bits,N) 
png(A_2,"a=2")

a=w0;
A_3=I_X_Vortex(w0,a,A,l,bits,N) 
png(A_3,"a=1")

a=.25*w0;
A_4=I_X_Vortex(w0,a,A,l,bits,N) 
png(A_4,"a=.25")

a=.5*w0;
A_5=I_X_Vortex(w0,a,A,l,bits,N) 
png(A_5,"a=.5")

#l variable, a fijo

l=.5;
L_1=I_X_Vortex(w0,a,A,l,bits,N)
png(L_1,"l=.5")

l=1;
L_2=I_X_Vortex(w0,a,A,l,bits,N)
png(L_2,"l=1")

l=2;
L_3=I_X_Vortex(w0,a,A,l,bits,N)
png(L_3,"l=2")

l=4;
L_4=I_X_Vortex(w0,a,A,l,bits,N)
png(L_4,"l=4")
