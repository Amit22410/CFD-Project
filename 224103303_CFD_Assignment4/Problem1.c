#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main ()
{
    int m1,n1;

    printf("Grid size m*n:\n m=\n n=\n  ");
    scanf("%d",&m1);
    scanf("%d",&n1);

    int m,n;
    m=m1+1;
    n=n1+1;
    double Re;

    printf("\nEnter Reynold's Number = ");
    scanf("%lf",&Re);

    double u[m][n+1], un[m][n+1], uc[m][n];
    double v[m+1][n], vn[m+1][n], vc[m][n];
    double p[m+1][n+1], pn[m+1][n+1], pc[m][n];
    double psi[m+2][n+2];
    double M[m+1][n+1];
    int i, j, step;
    double dx, dy, dt, delta, error;
    step =1;
    dx = 1.0/(m-1);
    dy = 1.0/(n-1);
    dt = 0.001;
    delta = 4.5;
    error = 1.0;

    for (i=0; i<=(m-1); i++)
    {

        for (j=0; j<=(n); j++)
        {

            u[i][j] = 0.0;
            u[i][n] = 1.0;
            u[i][n-1] = 1.0;
        }
    }

    for (i=0; i<=(m); i++)
    {

        for (j=0; j<=(n-1); j++)
        {

            v[i][j] = 0.0;
        }
    }

    for (i=0; i<=(m); i++)
    {

        for (j=0; j<=(n); j++)
        {

            p[i][j] = 1.0;
        }
    }

    while (error > 0.00001)
    {
        for (i=1; i<=(m-2); i++)
        {

            for (j=1; j<=(n-1); j++)
            {

                un[i][j] = u[i][j] - dt*(  (u[i+1][j]*u[i+1][j]-u[i-1][j]*u[i-1][j])/2.0/dx
                                           +0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i][j]+u[i][j-1])*(v[i+1][j-1]+v[i][j-1]) )/dy  )
                           - dt/dx*(p[i+1][j]-p[i][j])
                           + dt*1.0/Re*( (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx/dx +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy/dy );
            }
        }
        for (j=1; j<=(n-1); j++)
        {

            un[0][j] = 0.0;
            un[m-1][j] = 0.0;
        }

        for (i=0; i<=(n-1); i++)
        {

            un[i][0] = -un[i][1];
            un[i][n] = 2 - un[i][n-1];
        }

        for (i=1; i<=(m-1); i++)
        {

            for (j=1; j<=(n-2); j++)
            {

                vn[i][j] = v[i][j] - dt* ( 0.25*( (u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])-(u[i-1][j]+u[i-1][j+1])*(v[i][j]+v[i-1][j]) )/dx
                                           +(v[i][j+1]*v[i][j+1]-v[i][j-1]*v[i][j-1])/2.0/dy )
                           - dt/dy*(p[i][j+1]-p[i][j])
                           + dt*1.0/Re*( (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy/dy );
            }
        }
        for (j=1; j<=(n-2); j++)
        {

            vn[0][j] = -vn[1][j];
            vn[m][j] = -vn[m-1][j];
        }

        for (i=0; i<=(m); i++)
        {

            vn[i][0] = 0.0;
            vn[i][n-1] = 0.0;
        }

        for (i=1; i<=(m-1); i++)
        {

            for (j=1; j<=(n-1); j++)
            {

                pn[i][j] = p[i][j]-dt*delta*(  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] ) /dy  );
            }
        }

        for (i=1; i<=(m-1); i++)
        {

            pn[i][0] = pn[i][1];
            pn[i][n] = pn[i][n-1];
        }



        for (j=0; j<=(n); j++)
        {

            pn[0][j] = pn[1][j];
            pn[m][j] = pn[m-1][j];
        }


        error = 0.0;

        for (i=1; i<=(m-1); i++)
        {

            for (j=1; j<=(n-1); j++)
            {

                M[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );

                error = error + fabs(M[i][j]);
            }
        }

        if (step % 5000 ==1)
        {

            printf("\nError = %0.5lf for step = %d\n", error, step);

        }
        for (i=0; i<=(m-1); i++)
        {

            for (j=0; j<=(n); j++)
            {

                u[i][j] = un[i][j];
            }
        }
        for (i=0; i<=(m); i++)
        {

            for (j=0; j<=(n-1); j++)
            {

                v[i][j] = vn[i][j];
            }
        }

        for (i=0; i<=(m); i++)
        {

            for (j=0; j<=(n); j++)
            {

                p[i][j] = pn[i][j];
            }
        }
        step = step + 1;
    }

    for (i=0; i<=(m-1); i++)
    {

        for (j=0; j<=(n-1); j++)
        {

            uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
            vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
            pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
        }
    }

    FILE *ptr1, *ptr2, *ptr3, *ptr4, *ptr5 ;  //Output Files

    ptr1 = fopen("U_plot.dat","w+t");
    ptr2 = fopen("Central_U_plot.dat","w+t");
    ptr3 = fopen("V_plot.dat","w+t");
    ptr4 = fopen("Central_V_plot.dat","w+t");
    ptr5 = fopen("UV_streamline_plot.dat","w+t");

    fprintf( ptr1, "VARIABLES=\"X\",\"Y\",\"U\"\n");
    fprintf( ptr1, "ZONE  F=POINT\n");
    fprintf( ptr1, "I=%d, J=%d\n", m, n );

    for (j=0; j<(n); j++)
    {

        for (i=0; i<(m); i++ )
        {

            double xpos, ypos;
            xpos = i*dx;
            ypos = j*dy;

            fprintf( ptr1, "%.5lf\t%.5lf\t%.5lf\n", xpos, ypos, uc[i][j] );
        }
    }
    fprintf( ptr3, "VARIABLES=\"X\",\"Y\",\"V\"\n");
    fprintf( ptr3, "ZONE  F=POINT\n");
    fprintf( ptr3, "I=%d, J=%d\n", m, n );

    for (j=0; j<(n); j++)
    {

        for (i=0; i<(m); i++)
        {

            fprintf( ptr3, "%.5lf\t%.5lf\t%.5lf\n", i*dx, j*dy, vc[i][j] );
        }
    }
    fprintf(ptr2, "VARIABLES=\"U\",\"Y\"\n");
    fprintf(ptr2, "ZONE F=POINT\n");
    fprintf(ptr2, "J=%d\n", n );

    for (j=0; j<n; j++)
    {

        fprintf( ptr2, "%.5lf\t%.5lf\n", (uc[m/2][j] + uc[(m/2)+1][j])/(2.0),j*dy);
    }
    fprintf(ptr4, "VARIABLES=\"X\",\"V\"\n");
    fprintf(ptr4, "ZONE F=POINT\n");
    fprintf(ptr4, "I=%d\n", m);

    for (i=0; i<n; i++)
    {

        fprintf( ptr4, "%.5lf\t%.5lf\n",i*dx,(vc[i][n/2] + vc[i][n/2])/(2.) );
    }
    fprintf( ptr5, "VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n");
    fprintf( ptr5, "ZONE  F=POINT\n");
    fprintf( ptr5, "I=%d, J=%d\n", m, n );

    for (j=0; j<(n); j++)
    {

        for (i=0; i<(m); i++)
        {

            double xpos, ypos;
            xpos = i*dx;
            ypos = j*dy;

            fprintf( ptr5, "%.5lf\t%.5lf\t%.5lf\t%.5lf\n", xpos, ypos, uc[i][j], vc[i][j]);
        }
    }
    return 0;
}
