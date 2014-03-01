function printf(s,...)
    return io.write(s:format(...))
end

function fprintf(f,s,...)
    return f:write(s:format(...))
end

-- local versions of math functions
local sin, cos, sqrt, exp = math.sin, math.cos, math.sqrt, math.exp

--[[
acceleration functions and exact solutions
--]]

-- returns a function that only depends on position and velocity
function gravity()
    return function(vel, pos)
        return -9.81
    end
end

-- returns a function that only depends on absolute time
function gravityExact(acc0, vel0, pos0)
    return function(t)
        return vel0 + acc0 * t, pos0 + vel0 * t + acc0 * t * t * 0.5
    end
end

-- returns a function that only depends on position and velocity
function spring(stiffness, mass)
    return function(vel, pos)
        return -(stiffness / mass) * pos
    end
end

-- returns a function that only depends on absolute time
function springExact(stiffness, mass, acc0, vel0, pos0)
    local t0 = 0
    local denom = (sin(t0) ^ 2) + cos(t0)
    local A = (pos0 - vel0 * sin(t0)) / denom
    local B = (pos0 * sin(t0) + vel0 * cos(t0)) / denom
    -- print("initial conditions are A = ", A, "and B = ", B)
    return function(t)
        t = sqrt(stiffness / mass) * t
        local sint, cost = sin(t), cos(t)
        return -A * sint + B * cost, A * cost + B * sint
    end
end

-- returns a function that only depends on position and velocity
function dampedSpring(damping, stiffness, mass)
    return function(vel, pos)
        return -(damping * vel + stiffness * pos) / mass
    end
end

-- returns a function that only depends on absolute time
-- some constants that you would normally be able to choose are
-- fixed, because the exact solution is a mess...
function dampedSpringExact(damping, stiffness, mass, acc0, vel0, pos0)
    if vel0 ~= 0 or pos0 ~= 1 or mass ~= 1 then
        fprintf(io.stderr, "parameters are wrong: (vel0 = %f) != 0, (pos0 = %f) != 1, (mass = %f) != 1", vel0, pos0, mass)
        os.exit(1)
    end

    local c, k = damping, stiffness

    return function(t)
        -- t = sqrt(stiffness / mass) * t
        -- local sint, cost = sin(t), cos(t)
        local npos = (c + (c^2 - 4*k)^(1/2))/(2*exp(t*(c/2 - (c^2 - 4*k)^(1/2)/2))*(c^2 - 4*k)^(1/2)) -
            (c - (c^2 - 4*k)^(1/2))/(2*exp(t*(c/2 + (c^2 - 4*k)^(1/2)/2))*(c^2 - 4*k)^(1/2))
        return 0, npos
    end
end

--[[
Numerical integrators
--]]

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

do
    local old_acc = 0
    local have_old_acc = false
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
end

function ForwardEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)
    return vel + acc * dt, pos + vel * dt
end

-- some sort of midpoint method, as suggested by:
-- http://stackoverflow.com/questions/6446051/nsv-euler-integration-confusion?rq=1
-- http://gamedev.stackexchange.com/questions/25300/why-is-rk4-better-than-euler-integration
-- tries to keep the middle between symplectic Eulers
-- energy loss and forward Eulers energy addition
function MidpointEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)
    local hvel = vel + acc * dt * 0.5
    return vel + acc * dt, pos + hvel * dt
end

-- some sort of midpoint method, predictor-corrector style
-- should perform similar to real improved Euler, improves
-- the estimate of the velocity
function ImprovedMidpointEuler(accfn, vel, pos, dt)
    local acc = accfn(vel, pos)
    local hvel = vel + acc * dt * 0.5
    local hpos = vel + hvel * dt * 0.5
    acc = accfn(hvel, hpos)
    return vel + acc * dt, pos + hvel * dt
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

-- use a second-order method for approximating velocity
-- AND position. Should be the same as Velocity Verlet
-- according to my own checks. Vdrif says this is not
-- the case, let's see...
function ImprovedEulerVdrift(accfn, vel, pos, dt)
    local acc1 = accfn(vel, pos)
    local predpos = pos + vel * dt
    local predvel = vel + acc1 * dt

    local acc2 = accfn(predvel, predpos)
    local corrpos = pos + predvel * dt
    local corrvel = vel + acc2 * dt

    return corrvel, corrpos
end

function plotNumeric(timestep, it, integrator, accfn, pos, vel)
    for i = 1,it do
        printf("%d %f\n", i, pos)
        vel, pos = integrator(accfn, vel, pos, timestep(i))
    end
end

function plotExact(timestep, it, fn, accfn, pos0, vel0)
    for i = 1,it do
        local vel, pos = fn((i-1) * timestep(i))
        printf("%d %f\n", i, pos)
    end
end

function initPlot(name)
    printf("# X Y (%s)\n", name)
    print(name)
end

function nextPlot()
    io.write("\n\n")
end


local Gp = {}

function Gp:new(o)
    o = o or {}   -- create object if user does not provide one
    setmetatable(o, self)
    self.__index = self
    return o
end

function Gp:line(...)
    self.output:write(...)
    self.output:write("\n")
end

function Gp:title(title)
    self:line('set title "', title, '"')
end

function Gp:defineStyles()
    self:line("\n# define the styles")
    self:line("set style line 1 lc rgb '#0060ad' lw 2 pt 5 ps 1.5   # --- blue")
    self:line("set style line 2 lc rgb '#dd181f' lw 2 pt 7 ps 1.5   # --- red")
    self:line("set style line 3 lc rgb '#339900' lw 2 pt 9 ps 1.5   # --- green")
    self:line("set style line 4 lc rgb '#990066' lw 2 pt 11 ps 1.5   # --- magenta")
    self:line("set style line 5 lc rgb '#FF6633' lw 2 pt 13 ps 1.5   # --- orange")
    self:line("set style line 6 lc rgb '#262626' lw 2 pt 15 ps 1.5   # --- almost black")
    self:line("set style line 7 lc rgb '#599ad3' lw 2 pt 17 ps 1.5   # --- purple")
    self:line("set style line 8 lc rgb '#f9a65a' lw 2 pt 18 ps 1.5   # --- ???")
end

function Gp:plot(methods)
    local total = 0
    for _ in pairs(methods) do total = total + 1 end

    local str = "plot "
    local count = 0
    local template = "file index %d with linespoints ls %d"
    for name, fn in pairs(methods) do
        str = str .. template:format(count, count + 1)
        count = count + 1
        if total ~= count then
            str = str .. ", \\\n     "
        end
    end

    self:line("")
    self:line(str)
end

function Gp:init()
    self:line("set key autotitle columnhead")
    self:line("set grid")
end

function Gp:finish()
end

function plotMeta(out, title, methods)
    gp = Gp:new{output = out}
    gp:init()
    gp:title(title)
    gp:defineStyles()
    gp:plot(methods)
    gp:finish()
end

local methods = {}

-- first-order methods
methods["ForwardEuler"] = ForwardEuler
methods["SymplecticEuler"] = SymplecticEuler
methods["MidpointEuler"] = MidpointEuler

-- second-order methods
methods["RegularVerlet"] = Verlet
methods["TimeCorrectedVerlet"] = TimeCorrectedVerlet
methods["VelocityVerletVdrift"] = VelocityVerletVdrift
methods["NaiveImprovedEuler"] = NaiveImprovedEuler
methods["ImprovedEulerVdrift"] = ImprovedEulerVdrift

local compareToExact = false
local cmd = arg[1] or nil
if cmd == "data" then
    -- local timestep = function(iteration)
    --     return 1/60
    -- end
    -- local iterations = 20
    -- local vel0 = 1.6
    -- local pos0 = 5
    -- local acc = gravity()
    -- local exact = gravityExact(acc(vel0, pos0), vel0, pos0)

    -- local timestep = function(iteration)
    --     return 1/20
    -- end
    -- local iterations = 500
    -- local stiffness, mass = 1.5, 10
    -- local vel0, pos0 = 0, 1
    -- local acc = spring(stiffness, mass)
    -- local exact = springExact(stiffness, mass, acc(vel0, pos0), vel0, pos0)

    local timestep = function(iteration)
        return 1/30
    end
    local iterations = 500
    local drag, stiffness, mass = 0.5, 1.5, 1
    local vel0, pos0 = 0, 1
    local acc = dampedSpring(drag, stiffness, mass)
    local exact = dampedSpringExact(drag, stiffness, mass, acc(vel0, pos0), vel0, pos0)

    if compareToExact then
        -- plot the exact solution
       initPlot("Exact")
       plotExact(timestep, iterations, exact)
       nextPlot()
    end

    for name, fn in pairs(methods) do
        initPlot(name)
        plotNumeric(timestep, iterations, fn, acc, pos0, vel0)
        nextPlot()
    end
elseif cmd == "meta" then
    if compareToExact then
        methods["Exact"] = gravityExact
    end
    plotMeta(io.stdout, "numerical integration accuracy", methods)
else
    print("command ", cmd, " not recognized")
end

