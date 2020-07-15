function runge_kutta(func, x0, it, et, dt)
    
    t = range(it, et-dt, step=dt)
    b = et/dt
    b = convert(UInt64, b)
    time = reshape(t, (b,1))

    nt = size(t)
    nt = nt[1]

    nx = size(x0)
    nx = nx[1]

    x = zeros(nx,nt)

    for k in range(1, nt-1, step=1)
        k1 = dt .* func(t[k], x[:, k])
        k2 = dt .* func(t[k]+dt/2, x[:, k]+k1/2)
        k3 = dt .* func(t[k]+dt/2, x[:, k]+k2/2)
        k4 = dt .* func(t[k]+dt, x[:, k]+k3)

        test0 = k1
        test1 = 2*k2
        test2 = 2*k3
        test3 = k4

        dx = test0 + test1
        dx = dx + test2
        dx = dx + test3
        dx = dx ./ 6

        x[:, k+1] = x[:, k] + dx
    end

    return t,x
end