function printf(s,...)
    return io.write(s:format(...))
end

-- the acceleration function, best make this as complex as possible to
-- stress the integrators
function acceleration(vel, pos)
    return -9.81
end

function ConstantAccelerationExact(accfn, vel0, pos0, t)
    local acc0 = accfn()
    return vel0 + acc0 * t, pos0 + vel0 * t + acc0 * t * t * 0.5
end

old_acc = 0
have_old_acc = false
function VelocityVerletVdrift(accfn, vel, pos, dt)
    if have_old_acc ~= false then
        old_acc = accfn(vel, pos)
    end

    local npos = pos +  vel * dt + 0.5 * old_acc * dt * dt
    local nvel = vel +  0.5 * old_acc * dt

    -- calculate the acceleration at the end position and with the
    -- half-integrated velocity
    local acc = accfn(nvel, npos)

    -- correct the velocity
    nvel = nvel + 0.5 * acc * dt

    old_acc = acc
    have_old_acc = true

    return nvel, npos
end

function ForwardEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)
    return vel + acc * dt, pos + vel * dt
end

function SymplecticEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)

    -- forward euler step
    local nvel = vel + acc * dt

    -- backward euler step
    local npos = pos + nvel * dt

    return nvel, npos
end

function plotNumeric(it, integrator, accfn, pos, vel)
    for i = 1,it do
        printf("%d %f\n", i, pos)
        vel, pos = integrator(accfn, vel, pos, 1/60)
    end
end

function plotExact(it, fn, accfn, pos0, vel0)
    local vel = vel0
    local pos = pos0
    local t = 0
    for i = 1,it do
        printf("%d %f\n", i, pos)
        t = t + 1/60
        vel, pos = fn(accfn, vel0, pos0, t)
    end
end

function nextPlot()
    io.write("\n\n")
end

local iterations=20
local methods = {}

methods["Velocity Verlet Vdrift"] = VelocityVerletVdrift

function nextPlot()
    io.write("\n\n")
end

local iterations=20
local methods = {}

methods["VelocityVerletVdrift"] = VelocityVerletVdrift
methods["ForwardEuler"] = ForwardEuler
methods["SymplecticEuler"] = SymplecticEuler

printf("# X Y (exact)\n")
print("Exact")
plotExact(iterations, ConstantAccelerationExact, acceleration, 5, 1.6)
nextPlot()
for name, fn in pairs(methods) do
    printf("# X Y (%s)\n", name)
    print(name)
    plotNumeric(iterations, fn, acceleration, 5, 1.6)
    nextPlot()
end
