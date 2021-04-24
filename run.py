from methods import*


N = 50
vis_vel_field = True
dim = (N+2,N+2)
dt = 0.0004
diff = 0.01
visc = 0.01


u,v,u_prev,v_prev = [10*np.ones(dim) for i in range(4)]
dens,dens_prev,s = [np.zeros(dim) for i in range(3)]



su,sv = [0*np.ones(dim) for i in range(2)]
s += 0.1

su[:10,:10] += 10000
sv[:10,:10] += 7000

set_bnd(N,1,u)
set_bnd(N,2,v)

#su = -su
#sv = -sv

tmax = 100000*dt





size = int(np.ceil(15*50/N))
framerate = 60
frametime = 1/framerate
fontsize = 25









pygame.init()

gridDisplay = pygame.display.set_mode(((N+2)*size, (N+2)*size))
pygame.display.get_surface().fill((200, 200, 200))  # background


# we use the sizes to draw as well as to do our "steps" in the loops.
grid_node_width = size
grid_node_height = size


pygame.font.init() # you have to call this at the start,
                   # if you want to use this module.
myfont = pygame.font.SysFont('lato', fontsize,bold = True)

def createSquare(x, y, color):
    pygame.draw.rect(gridDisplay, color, [x, y, grid_node_width, grid_node_height ])

def createArrow(x,y,x1,y1):
    corr_x = x + 1/2*grid_node_width
    corr_y = y + 1/2*grid_node_height
    pygame.draw.line(gridDisplay, (255,255,255), (corr_x,corr_y), (corr_x+x1, corr_y+y1))
    pygame.draw.circle(gridDisplay, (255,255,255), (corr_x+x1,corr_y+y1),size/7)

def visualizeGrid(mat,u,v,vis_vel_field = True):
    maxi = mat.max()
    mini = mat.min()
    matrix = 255*(mat-mini)/(maxi-mini)

    ubermax = max(max(abs(u.max()),abs(u.min())),max(abs(v.max()),abs(v.min())))
    u = 1.5*1/2*grid_node_height*u/ubermax
    v = 1.5*1/2*grid_node_width*v/ubermax
    y = 0  # we start at the top of the screen
    for i in range(N+2):
        x = 0 # for every row we start at the left of the screen again
        for j in range(N+2):

            createSquare(x, y, (matrix[i,j],0,0))
            x += grid_node_width # for ever item/number in that row we move one "step" to the right
        y += grid_node_height   # for every new row we move one "step" downwards
    if vis_vel_field:
        y = 0  # we start at the top of the screen
        for i in range(N+2):
            x = 0 # for every row we start at the left of the screen again
            for j in range(N+2):
                createArrow(x,y,v[i,j],u[i,j])
                x += grid_node_width # for ever item/number in that row we move one "step" to the right
            y += grid_node_height   # for every new row we move one "step" downwards

    return ubermax



t = 0
prev_start = 0
while t < tmax:
    start = time.time()
    #if t > 5*dt:
        #s = 0*s
    if t > 600*dt:
        su *= 0
        sv *= 0
    vel_step(N,u,v,u_prev,v_prev,visc,dt,su,sv)
    dens_step(N,dens,dens_prev,u,v,diff,dt,s)
    ubermax = visualizeGrid(dens,u,v,vis_vel_field)
    textsurface = myfont.render('FPS: ' + str(np.round(1/(start-prev_start),5)), False, (0, 255, 50))
    textsurface1 = myfont.render('Velocity: ' + str(np.round(ubermax,5)), False, (0, 255, 50))
    textsurface2 = myfont.render('Concentration: ' + str(np.round(np.sum(dens),5)), False, (0, 255, 50))
    gridDisplay.blit(textsurface,(0,0))
    gridDisplay.blit(textsurface1,(0,0.035*(N+2)*size*fontsize/33))
    gridDisplay.blit(textsurface2,(0,0.07*(N+2)*size*fontsize/33))
    pygame.display.update()
    while time.time() -start < frametime:
        pass
    t+=dt
    prev_start = start


 # call the function

start = time.time()

while time.time() - start < 1:
    pass

pygame.font.quit()
pygame.quit()
#quit()
