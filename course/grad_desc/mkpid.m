function C = mkpid(th)
C = pid(exp(th(1)), exp(th(2)), exp(th(3)), exp(th(4)));
end