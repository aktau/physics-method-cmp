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

-- use a second-order method for approximating velocity
-- AND position. Should be the same as Velocity Verlet
-- according to my own checks. Vdrif says this is not
-- the case, let's see...
function ImprovedEuler(accfn, vel, pos, dt)
    local acc1 = accfn(vel, pos)
    local predpos = pos + vel * dt
    local predvel = vel + acc1 * dt

    local acc2 = accfn(predvel, predpos)
    local corrpos = pos + predvel * dt
    local corrvel = vel + acc2 * dt

    return corrvel, corrpos
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


local Gp = {}

function Gp:new(o)
    o = o or {}   -- create object if user does not provide one
    setmetatable(o, self)
    self.__index = self
    return o
end

function Gp:line(...)
    -- print(...)
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
    self:line("set style line 7 lc rgb '#599ad3' lw 2 pt 15 ps 1.5   # --- purple")
    self:line("set style line 8 lc rgb '#f9a65a' lw 2 pt 15 ps 1.5   # --- ???")
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

local iterations=20
local methods = {}

methods["VelocityVerletVdrift"] = VelocityVerletVdrift
methods["TimeCorrectedVerlet"] = TimeCorrectedVerlet
methods["ForwardEuler"] = ForwardEuler
methods["SymplecticEuler"] = SymplecticEuler
methods["RegularVerlet"] = Verlet
methods["NaiveImprovedEuler"] = NaiveImprovedEuler
methods["ImprovedEuler"] = ImprovedEuler

local cmd = arg[1] or nil
if cmd == "data" then
    -- plot the exact solution
    initPlot("Exact")
    plotExact(iterations, ConstantAccelerationExact, acceleration, 5, 1.6)
    nextPlot()

    for name, fn in pairs(methods) do
        initPlot(name)
        plotNumeric(iterations, fn, acceleration, 5, 1.6)
        nextPlot()
    end
elseif cmd == "meta" then
    plotMeta(io.stdout, "numerical integration accuracy", methods)
else
    print("command ", cmd, " not recognized")
end

