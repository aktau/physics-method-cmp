function printf(s,...)
    return io.write(s:format(...))
end

-- the acceleration function, best make this as complex as possible to
-- stress the integrators
function acceleration(vel, pos)
    return -9.81
end

old_acc = 0
have_old_acc = false
function VelocityVerletVdrift(vel, pos, dt)
    if have_old_acc ~= false then
        old_acc = acceleration(vel, pos)
    end

    local npos = pos +  vel * dt + 0.5 * old_acc * dt * dt
    local nvel = vel +  0.5 * old_acc * dt

    -- calculate the acceleration at the end position and with the
    -- half-integrated velocity
    local acc = acceleration(nvel, npos)

    -- correct the velocity
    nvel = nvel + 0.5 * acc * dt

    old_acc = acc
    have_old_acc = true

    return nvel, npos
end

function ForwardEuler(vel, pos, dt)
    local acc = acceleration(vel, pos)
    return vel + acc * dt, pos + vel * dt
end

function plot(fn, it, pos, vel)
    print("# X Y")
    for i = 1,it do
        -- print(pos, vel)
        printf("%d %f\n", i, pos)
        pos, vel = fn(vel, pos, 1/60)
    end
end

function nextPlot()
    io.write("\n\n")
end

plot(VelocityVerletVdrift, 5, 5, 1.6)
nextPlot()
plot(ForwardEuler, 5, 5, 1.6)
