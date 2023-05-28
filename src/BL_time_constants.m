function [Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,qR,R,RD,P,qR_max,lambda_1,lambda_2,Ta,Tf)

%% Stall supression time delay 
Ta_theta = Ta*(1/4+lambda_1*RD^(lambda_2)*(1-R^(1/4))*exp(-(abs(theta)-1)^2/0.01));

%% Separation points time delays         
Tf_n = Tf*(1+1/2*P^6*~upstroke+RD^3*(4*(1-qR/qR_max)*~upstroke+1/2*(qR/qR_max)^4*upstroke));    % +1/2*P^6*~upstroke to increase the delay on the downstroke of very light stall, +RD^3*(4*(1-qR/qR_max)*~upstroke+1/2*(qR/qR_max)^4*upstroke to increase the delay at the end of the downstroke and begin of upstroke for very high pitch rates
Tf_m = Tf*(1-1/2*R*(1-P)*~upstroke);                                                            % -1/2*R*(1-P)*~upstroke to reduce delay at high pitch rates in the downstroke, but not in light stall
Tf_c = Tf_n;  

end