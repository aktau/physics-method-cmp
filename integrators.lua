#!/usr/bin/env luajit

printf = function(s,...)
    return io.write(s:format(...))
end -- function

function acceleration(vel, pos)
    return -9.81
end

-- if (!have_oldaccel)
--     oldaccel = system.GetAcceleration(state);

-- state.x += state.v*dt + 0.5*oldaccel*dt*dt;
-- state.v += 0.5*oldaccel*dt;
-- real a = system.GetAcceleration(state);
-- state.v += 0.5*a*dt;

-- oldaccel = a;
-- have_oldaccel = true;
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

function plot(fn, it, pos, vel)
    print("# X Y")
    for i = 1,it do
        -- print(pos, vel)
        printf("%d %f\n", i, pos)
        pos, vel = fn(pos, vel, 1/60)
    end
end

plot(VelocityVerletVdrift, 5, 5, 1.6)
