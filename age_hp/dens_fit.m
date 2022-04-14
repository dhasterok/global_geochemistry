test = data(data.density > 0 & data.density < 4000 & ...
    data.sio2 >= 0 & data.sio2 <= 100 &...
    data.al2o3 >= 0 & data.al2o3 <= 100 &...
    data.feo_tot >= 0 & data.feo_tot <= 100 &...
    data.mgo >= 0 & data.mgo <= 100 &...
    data.cao >= 0 & data.cao <= 100 &...
    data.na2o >= 0 & data.na2o <= 100 &...
    data.tio2 >= 0 & data.tio2 <= 100 &...
    data.k2o >= 0 & data.k2o <= 100 &...
    data.p2o5 >= 0 & data.p2o5 <= 100,:);

ind = find(test.density < 10);
test.density(ind) = test.density(ind) * 1000;




y = test.density;

fun = @(x)x(1).*abs(test.sio2) + x(2).*abs(test.al2o3) + x(3).*abs(test.feo_tot)...
    + x(4).*abs(test.mgo) + x(4).*abs(test.cao) - y;

sio2_guess = 5;
al2o3_guess = 3;
feo_tot_guess = 1;
mgo_guess = 1;
cao_guess = 1;
na2o_guess = 1;
tio2_guess = 1;
k2o_guess = 1;
p2o5_guess = 1;

start_point(1,1) = sio2_guess;
start_point(1,2) = al2o3_guess;
start_point(1,3) = feo_tot_guess;
start_point(1,4) = mgo_guess;
start_point(1,5) = cao_guess;
start_point(1,6) = na2o_guess;
start_point(1,7) = tio2_guess;
start_point(1,8) = k2o_guess;
start_point(1,9) = p2o5_guess;
st=start_point;
x0 = st; % arbitrary initial guess

x = lsqnonlin(fun,x0);

plot(y,fun(x)+y,'r.')
hold on
plot([2000:200:4000],[2000:200:4000],'-g')
legend('Data','Best fit')
xlabel('Measured density')
ylabel('Estimated density')
xlim([2400,3400])
ylim([2400,3400])
hold off