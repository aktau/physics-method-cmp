function printf(s,...)
    return io.write(s:format(...))
end

function fprintf(f,s,...)
    return f:write(s:format(...))
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

do
    local old_pos = 0
    local have_old_pos = false

    function Verlet(accfn, vel, pos, dt)
        if have_old_pos == false then
            -- first pass is not true Verlet
            have_old_pos = true

            local npos = pos + vel * dt + accfn(vel, pos) * dt * dt * 0.5
            old_pos = pos
            return 0, npos
        end

        local npos = pos + (pos - old_pos) + accfn(vel, pos) * dt * dt
        old_pos = pos
        return 0, npos
    end
end

do
    local old_pos = 0
    local old_dt = 0
    local have_old_pos = false

    -- based on
    -- http://archive.gamedev.net/archive/reference/programming/features/verlet/default.html
    function TimeCorrectedVerlet(accfn, vel, pos, dt)
        if have_old_pos == false then
            -- first pass is not true Verlet
            have_old_pos = true

            local npos = pos + vel * dt + accfn(vel, pos) * dt * dt * 0.5
            old_pos = pos
            old_dt = dt
            return 0, npos
        end

        local npos = pos + (pos - old_pos) * (dt / old_dt) + accfn(vel, pos) * dt * dt
        old_pos = pos
        old_dt = dt
        return 0, npos
    end
end

old_acc = 0
have_old_acc = false
function VelocityVerletVdrift(accfn, vel, pos, dt)
    if have_old_acc == false then
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
    local nvel = vel + acc * dt -- forward euler step
    local npos = pos + nvel * dt -- backward euler step
    return nvel, npos
end

-- use a second-order method for approximating position
-- should give the same accuracy as the Verlet family for
-- when using constant acceleration
function NaiveImprovedEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)
    local nvel = vel + acc * dt
    return nvel, pos + (vel + nvel) * 0.5 * dt
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

function initPlot(name)
    printf("# X Y (%s)\n", name)
    print(name)
end

function nextPlot()
    io.write("\n\n")
end

local iterations=20
local methods = {}

methods["VelocityVerletVdrift"] = VelocityVerletVdrift
methods["TimeCorrectedVerlet"] = TimeCorrectedVerlet
methods["ForwardEuler"] = ForwardEuler
methods["SymplecticEuler"] = SymplecticEuler
methods["RegularVerlet"] = Verlet
methods["NaiveImprovedEuler"] = NaiveImprovedEuler

initPlot("Exact")
plotExact(iterations, ConstantAccelerationExact, acceleration, 5, 1.6)
nextPlot()
for name, fn in pairs(methods) do
    initPlot(name)
    plotNumeric(iterations, fn, acceleration, 5, 1.6)
    nextPlot()
end
